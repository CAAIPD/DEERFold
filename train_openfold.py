import argparse
import logging
import os, sys
import gc

#os.environ["CUDA_VISIBLE_DEVICES"] = "0"
#os.environ["MASTER_ADDR"]="10.119.81.14"
#os.environ["MASTER_PORT"]="42069"
#os.environ["NODE_RANK"]="0"

import random
import time

import numpy as np
import pytorch_lightning as pl
from pytorch_lightning.callbacks.lr_monitor import LearningRateMonitor
from pytorch_lightning.callbacks.model_checkpoint import ModelCheckpoint
from pytorch_lightning.plugins.training_type import DeepSpeedPlugin, DDPPlugin
from pytorch_lightning.plugins.environments import SLURMEnvironment
from pytorch_lightning.loggers import WandbLogger
import torch

from openfold.config import model_config
from openfold.data.data_modules import (
    OpenFoldDataModule,
    DummyDataLoader,
)
from openfold.model.model import AlphaFold
from openfold.model.torchscript import script_preset_
from openfold.np import residue_constants
from openfold.utils.callbacks import (
    EarlyStoppingVerbose,
)
from openfold.utils.exponential_moving_average import ExponentialMovingAverage
from openfold.utils.argparse import remove_arguments
from openfold.utils.loss import AlphaFoldLoss, lddt_ca, compute_drmsd
from openfold.utils.seed import seed_everything
from openfold.utils.superimposition import superimpose
from openfold.utils.tensor_utils import tensor_tree_map
from openfold.utils.validation_metrics import (
    drmsd,
    gdt_ts,
    gdt_ha,
)
from scripts.zero_to_fp32 import (
    get_fp32_state_dict_from_zero_checkpoint
)

from openfold.utils.import_weights import (
    import_jax_weights_,
)

from openfold.utils.logger import PerformanceLoggingCallback

def filter_batch(batch):
    keys_to_keep = [
        'aatype', 'residue_index', 'seq_length', 'sp', 'seq_mask', 'msa_mask',
        'msa_row_mask', 'atom14_atom_exists', 'residx_atom14_to_atom37',
        'residx_atom37_to_atom14', 'atom37_atom_exists', 'extra_msa',
        'extra_msa_mask', 'extra_msa_row_mask', 'bert_mask', 'true_msa',
        'extra_has_deletion', 'extra_deletion_value', 'msa_feat', 'target_feat', 'all_atom_positions', 'all_atom_mask'
    ]
    
    return {key: value for key, value in batch.items() if key in keys_to_keep}

class OpenFoldWrapper(pl.LightningModule):
    def __init__(self, config, model):
        super(OpenFoldWrapper, self).__init__()
        self.config = config
        self.model = model #AlphaFold(config)
        self.loss = AlphaFoldLoss(config.loss)
        self.ema = ExponentialMovingAverage(
            model=self.model, decay=config.ema.decay
        )
        
        self.cached_weights = None

    def forward(self, batch):
        return self.model(batch)

    def _log(self, loss_breakdown, batch, outputs, train=True):
        phase = "train" if train else "val"
        for loss_name, indiv_loss in loss_breakdown.items():
            self.log(
                f"{phase}/{loss_name}", 
                indiv_loss, 
                on_step=train, on_epoch=(not train), logger=True,
            )

            if(train):
                self.log(
                    f"{phase}/{loss_name}_epoch",
                    indiv_loss,
                    on_step=False, on_epoch=True, logger=True,
                )

        with torch.no_grad():
            other_metrics = self._compute_validation_metrics(
                batch, 
                outputs,
                superimposition_metrics=(not train)
            )

        for k,v in other_metrics.items():
            self.log(
                f"{phase}/{k}",
                torch.mean(v),
                on_step=False, on_epoch=True, logger=True
            )

    def training_step(self, batch, batch_idx):
        if(self.ema.device != batch["aatype"].device):
            self.ema.to(batch["aatype"].device)

        # Run the model
        outputs = self(batch)
        
        # Remove the recycling dimension
        batch = tensor_tree_map(lambda t: t[..., -1], batch)

        # Compute loss
        loss, loss_breakdown = self.loss(
            outputs, batch, _return_breakdown=True
        )

        # Log it
        self._log(loss_breakdown, batch, outputs)

        return loss

    def on_validation_start(self):
        self.model.eval()
    
    def validation_step(self, batch, batch_idx):
        # At the start of validation, load the EMA weights
        if(self.cached_weights is None):
            self.cached_weights = self.model.state_dict()
            self.model.load_state_dict(self.ema.state_dict()["params"])
        
        # Calculate validation loss
        batch = filter_batch(batch)

        try:
                outputs = self.model(batch)
                batch = tensor_tree_map(lambda t: t[..., -1], batch)
                metrics = self._compute_validation_metrics(
                    batch, 
                    outputs,
                    superimposition_metrics=True
                )
                for k, v in metrics.items():
                    # print(f"val/{k}", v)
                    self.log(f"val/{k}", v)

                # loss = lddt_ca(
                #     outputs["final_atom_positions"],
                #     batch["all_atom_positions"],
                #     batch["all_atom_mask"],
                #     eps=self.config.globals.eps,
                #     per_residue=False,
                # )

                self.log("val_loss", metrics['gdt_ts'], prog_bar=True)
                del outputs, batch
                torch.cuda.empty_cache()
                gc.collect()
                return {"val_loss": metrics['gdt_ts']}
            
        except torch.cuda.OutOfMemoryError:
            print(f"OOM error occurred for batch {batch_idx}. Skipping this step.")
            torch.cuda.empty_cache() 
            return {"val_loss": None}
        
        except Exception as e:
            print(f"An error occurred for batch {batch_idx}: {str(e)}. Skipping this step.")
            return {"val_loss": None} 

    def validation_epoch_end(self, _):
        # Restore the model weights to normal
        self.model.load_state_dict(self.cached_weights)
        self.cached_weights = None

    def _compute_validation_metrics(self, 
        batch, 
        outputs, 
        superimposition_metrics=False
    ):
        metrics = {}
        
        gt_coords = batch["all_atom_positions"]
        pred_coords = outputs["final_atom_positions"]
        all_atom_mask = batch["all_atom_mask"]
    
        # This is super janky for superimposition. Fix later
        gt_coords_masked = gt_coords * all_atom_mask[..., None]
        pred_coords_masked = pred_coords * all_atom_mask[..., None]
        ca_pos = residue_constants.atom_order["CA"]
        gt_coords_masked_ca = gt_coords_masked[..., ca_pos, :]
        pred_coords_masked_ca = pred_coords_masked[..., ca_pos, :]
        all_atom_mask_ca = all_atom_mask[..., ca_pos]
    
        lddt_ca_score = lddt_ca(
            pred_coords,
            gt_coords,
            all_atom_mask,
            eps=self.config.globals.eps,
            per_residue=False,
        )
        # print(batch['seq_length'], lddt_ca_score)
   
        metrics["lddt_ca"] = lddt_ca_score

        # drmsd_ca_score = drmsd(
        #     pred_coords_masked_ca,
        #     gt_coords_masked_ca,
        #     mask=all_atom_mask_ca, # still required here to compute n
        # )
   
        # metrics["drmsd_ca"] = drmsd_ca_score
    
        if(superimposition_metrics):
            superimposed_pred, alignment_rmsd = superimpose(
                gt_coords_masked_ca, pred_coords_masked_ca, all_atom_mask_ca,
            )
            gdt_ts_score = gdt_ts(
                superimposed_pred, gt_coords_masked_ca, all_atom_mask_ca
            )
            gdt_ha_score = gdt_ha(
                superimposed_pred, gt_coords_masked_ca, all_atom_mask_ca
            )

            metrics["alignment_rmsd"] = alignment_rmsd
            metrics["gdt_ts"] = gdt_ts_score
            metrics["gdt_ha"] = gdt_ha_score
    
        return metrics

    def configure_optimizers(self, 
        learning_rate: float = 1e-3,
        eps: float = 1e-8
    ) -> torch.optim.Adam:
        # Ignored as long as a DeepSpeed optimizer is configured
        return torch.optim.Adam(
            self.model.parameters(), 
            lr=learning_rate, 
            eps=eps
        )

    def on_before_zero_grad(self, *args, **kwargs):
        self.ema.update(self.model)

    def on_save_checkpoint(self, checkpoint):
        checkpoint["ema"] = self.ema.state_dict()


def main(args):
    if(args.seed is not None):
        seed_everything(args.seed) 

    config = model_config(
        "model_5_ptm", 
        train=True, 
        low_prec=(args.precision == 16)
    ) 
    

    model = AlphaFold(config)   
    
    import_jax_weights_(model, "/data/mchaourab/wut18/release/DEERFold/openfold/resources/params/params_model_5_ptm.npz", version="model_5_ptm")

    model_module = OpenFoldWrapper(config, model)

    #script_preset_(model)

    if(args.resume_from_ckpt and args.resume_model_weights_only):
        sd = get_fp32_state_dict_from_zero_checkpoint(args.resume_from_ckpt)
        sd = {k[len("module."):]:v for k,v in sd.items()}
        model_module.load_state_dict(sd)
        logging.info("Successfully loaded model weights...")

    # TorchScript components of the model
    if(args.script_modules):
        script_preset_(model_module)

    #data_module = DummyDataLoader("batch.pickle")
    data_module = OpenFoldDataModule(
        config=config.data, 
        batch_seed=args.seed,
        **vars(args)
    )

    data_module.prepare_data()
    data_module.setup()
    
    callbacks = []
    mc = ModelCheckpoint(
        dirpath=os.path.join(args.output_dir, "checkpoints"),
        filename="deerfold_{epoch}_{val_loss:.4f}",
        every_n_epochs=1,
        save_top_k=-1,
    )
    callbacks.append(mc)

    if(args.early_stopping):
        es = EarlyStoppingVerbose(
            monitor="val_loss",
            min_delta=args.min_delta,
            patience=args.patience,
            verbose=False,
            mode="max",
            check_finite=True,
            strict=True,
        )
        callbacks.append(es)
        
    if(args.log_performance):
        global_batch_size = args.num_nodes * args.gpus
        perf = PerformanceLoggingCallback(
            log_file=os.path.join(args.output_dir, "performance_log.json"),
            global_batch_size=global_batch_size,
        )
        callbacks.append(perf)

    if(args.log_lr):
        lr_monitor = LearningRateMonitor(logging_interval="step")
        callbacks.append(lr_monitor)

    loggers = []
    if(args.wandb):
        wdb_logger = WandbLogger(
            name=args.experiment_name,
            save_dir=args.output_dir,
            id=args.wandb_id,
            project=args.wandb_project,
            **{"entity": args.wandb_entity}
        )
        loggers.append(wdb_logger)

    if(args.deepspeed_config_path is not None):
        strategy = DeepSpeedPlugin(
            config=args.deepspeed_config_path,
        )
        if(args.wandb):
            wdb_logger.experiment.save(args.deepspeed_config_path)
            wdb_logger.experiment.save("openfold/config.py")
    elif (args.gpus is not None and args.gpus > 1) or args.num_nodes > 1:
        strategy = DDPPlugin(find_unused_parameters=False)
    else:
        strategy = None
 
    if(args.wandb):
        freeze_path = f"{wdb_logger.experiment.dir}/package_versions.txt"
        os.system(f"{sys.executable} -m pip freeze > {freeze_path}")
        wdb_logger.experiment.save(f"{freeze_path}")

    trainer = pl.Trainer.from_argparse_args(
        args,
        default_root_dir=args.output_dir,
        strategy=strategy,
        callbacks=callbacks,
        logger=loggers,
    )

    if(args.resume_model_weights_only):
        ckpt_path = None
    else:
        ckpt_path = args.resume_from_ckpt

    trainer.fit(
        model_module, 
        datamodule=data_module,
        ckpt_path=ckpt_path,
    )    

    save_dir =  wdb_logger.save_dir
    checkpoint_path = os.path.join(save_dir, wdb_logger.name, "checkpoints", "final.ckpt")
    print(checkpoint_path)
    trainer.save_checkpoint(checkpoint_path)
    # trainer.save_checkpoint(
    #     os.path.join(trainer.logger.log_dir, "checkpoints", "final.ckpt")
    # )

def bool_type(bool_str: str):
    bool_str_lower = bool_str.lower()
    if bool_str_lower in ('false', 'f', 'no', 'n', '0'):
        return False
    elif bool_str_lower in ('true', 't', 'yes', 'y', '1'):
        return True
    else:
        raise ValueError(f'Cannot interpret {bool_str} as bool')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "train_data_dir", type=str,
        help="Directory containing training mmCIF files"
    )
    parser.add_argument(
        "train_alignment_dir", type=str,
        help="Directory containing precomputed training alignments"
    )
    parser.add_argument(
        "template_mmcif_dir", type=str,
        help="Directory containing mmCIF files to search for templates"
    )
    parser.add_argument(
        "output_dir", type=str,
        help='''Directory in which to output checkpoints, logs, etc. Ignored
                if not on rank 0'''
    )
    parser.add_argument(
        "--max_template_date", type=str,  default="2021-10-10",
        help='''Cutoff for all templates. In training mode, templates are also 
                filtered by the release date of the target'''
    )
    parser.add_argument(
        "--train_feats_dir", type=str, default=None,
        help="Directory containing training feature files in pkl"
    )
    parser.add_argument(
        "--train_DEER_dir", type=str,
        help="Directory containing training DEER data in csv"
    )
    parser.add_argument(
        "--distillation_data_dir", type=str, default=None,
        help="Directory containing training PDB files"
    )
    parser.add_argument(
        "--distillation_alignment_dir", type=str, default=None,
        help="Directory containing precomputed distillation alignments"
    )
    parser.add_argument(
        "--val_data_dir", type=str, default=None,
        help="Directory containing validation mmCIF files"
    )
    parser.add_argument(
        "--val_feats_dir", type=str, default=None,
        help="Directory containing validation feature files in pkl"
    )
    parser.add_argument(
        "--val_DEER_dir", type=str,
        help="Directory containing validation DEER data in csv"
    )
    parser.add_argument(
        "--val_alignment_dir", type=str, default=None,
        help="Directory containing precomputed validation alignments"
    )
    parser.add_argument(
        "--kalign_binary_path", type=str, default='/usr/bin/kalign',
        help="Path to the kalign binary"
    )
    parser.add_argument(
        "--train_mapping_path", type=str, default=None,
        help='''Optional path to a .json file containing a mapping from
                consecutive numerical indices to sample names. Used to filter
                the training set'''
    )
    parser.add_argument(
        "--val_mapping_path", type=str, default=None,
        help='''Optional path to a .json file containing a mapping from
                consecutive numerical indices to sample names. Used to filter
                the validation set'''
    )
    parser.add_argument(
        "--distillation_mapping_path", type=str, default=None,
        help="""See --train_mapping_path"""
    )
    parser.add_argument(
        "--obsolete_pdbs_file_path", type=str, default=None,
        help="""Path to obsolete.dat file containing list of obsolete PDBs and 
             their replacements."""
    )
    parser.add_argument(
        "--template_release_dates_cache_path", type=str, default=None,
        help="""Output of scripts/generate_mmcif_cache.py run on template mmCIF
                files."""
    )
    parser.add_argument(
        "--use_small_bfd", type=bool_type, default=False,
        help="Whether to use a reduced version of the BFD database"
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed"
    )
    parser.add_argument(
        "--deepspeed_config_path", type=str, default=None,
        help="Path to DeepSpeed config. If not provided, DeepSpeed is disabled"
    )
    parser.add_argument(
        "--checkpoint_best_val", type=bool_type, default=True,
        help="""Whether to save the model parameters that perform best during
                validation"""
    )
    parser.add_argument(
        "--early_stopping", type=bool_type, default=False,
        help="Whether to stop training when validation loss fails to decrease"
    )
    parser.add_argument(
        "--min_delta", type=float, default=0,
        help="""The smallest decrease in validation loss that counts as an 
                improvement for the purposes of early stopping"""
    )
    parser.add_argument(
        "--patience", type=int, default=3,
        help="Early stopping patience"
    )
    parser.add_argument(
        "--resume_from_ckpt", type=str, default=None,
        help="Path to a model checkpoint from which to restore training state"
    )
    parser.add_argument(
        "--resume_model_weights_only", type=bool_type, default=False,
        help="Whether to load just model weights as opposed to training state"
    )
    parser.add_argument(
        "--log_performance", type=bool_type, default=False,
        help="Measure performance"
    )
    parser.add_argument(
        "--log_lr", action="store_true", default=False,
        help="Whether to log the actual learning rate"
    )
    parser.add_argument(
        "--wandb", action="store_true", default=False,
        help="Whether to log metrics to Weights & Biases"
    )
    parser.add_argument(
        "--script_modules", type=bool_type, default=False,
        help="Whether to TorchScript eligible components of them model"
    )
    parser.add_argument(
        "--experiment_name", type=str, default=None,
        help="Name of the current experiment. Used for wandb logging"
    )
    parser.add_argument(
        "--wandb_id", type=str, default=None,
        help="ID of a previous run to be resumed"
    )
    parser.add_argument(
        "--wandb_project", type=str, default=None,
        help="Name of the wandb project to which this run will belong"
    )
    parser.add_argument(
        "--wandb_entity", type=str, default=None,
        help="wandb username or team name to which runs are attributed"
    )

    parser = pl.Trainer.add_argparse_args(parser)
   
    # Disable the initial validation pass
    parser.set_defaults(
        num_sanity_val_steps=0,
    )

    # Remove some buggy/redundant arguments introduced by the Trainer
    remove_arguments(parser, ["--accelerator", "--resume_from_checkpoint"]) 

    args = parser.parse_args()

    if(args.seed is None and 
        ((args.gpus is not None and args.gpus > 1) or 
         (args.num_nodes is not None and args.num_nodes > 1))):
        raise ValueError("For distributed training, --seed must be specified")

    main(args)
