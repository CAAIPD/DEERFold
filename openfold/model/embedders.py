# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import torch
import torch.nn as nn
from typing import Tuple

from openfold.model.primitives import Linear
from openfold.utils.tensor_utils import one_hot


class InputEmbedder(nn.Module):
    """
    Embeds a subset of the input features.

    Implements Algorithms 3 (InputEmbedder) and 4 (relpos).
    """

    def __init__(
        self,
        tf_dim: int,
        msa_dim: int,
        c_z: int,
        c_m: int,
        relpos_k: int,
        **kwargs,
    ):
        """
        Args:
            tf_dim:
                Final dimension of the target features
            msa_dim:
                Final dimension of the MSA features
            c_z:
                Pair embedding dimension
            c_m:
                MSA embedding dimension
            relpos_k:
                Window size used in relative positional encoding
        """
        super(InputEmbedder, self).__init__()

        self.tf_dim = tf_dim
        self.msa_dim = msa_dim

        self.c_z = c_z
        self.c_m = c_m

        self.linear_tf_z_i = Linear(tf_dim, c_z)
        self.linear_tf_z_j = Linear(tf_dim, c_z)
        self.linear_tf_m = Linear(tf_dim, c_m)
        self.linear_msa_m = Linear(msa_dim, c_m)

        # RPE stuff
        self.relpos_k = relpos_k
        self.no_bins = 2 * relpos_k + 1
        self.linear_relpos = Linear(self.no_bins, c_z)

    def relpos(self, ri: torch.Tensor):
        """
        Computes relative positional encodings

        Implements Algorithm 4.

        Args:
            ri:
                "residue_index" features of shape [*, N]
        """
        d = ri[..., None] - ri[..., None, :]
        boundaries = torch.arange(
            start=-self.relpos_k, end=self.relpos_k + 1, device=d.device
        )
        oh = one_hot(d, boundaries).type(ri.dtype)
        return self.linear_relpos(oh)

    def forward(
        self,
        tf: torch.Tensor,
        ri: torch.Tensor,
        msa: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Args:
            tf:
                "target_feat" features of shape [*, N_res, tf_dim]
            ri:
                "residue_index" features of shape [*, N_res]
            msa:
                "msa_feat" features of shape [*, N_clust, N_res, msa_dim]
        Returns:
            msa_emb:
                [*, N_clust, N_res, C_m] MSA embedding
            pair_emb:
                [*, N_res, N_res, C_z] pair embedding

        """
        # [*, N_res, c_z]
        tf_emb_i = self.linear_tf_z_i(tf)
        tf_emb_j = self.linear_tf_z_j(tf)

        # [*, N_res, N_res, c_z]
        pair_emb = tf_emb_i[..., None, :] + tf_emb_j[..., None, :, :]
        pair_emb = pair_emb + self.relpos(ri.type(pair_emb.dtype))

        # [*, N_clust, N_res, c_m]
        n_clust = msa.shape[-3]
        tf_m = (
            self.linear_tf_m(tf)
            .unsqueeze(-3)
            .expand(((-1,) * len(tf.shape[:-2]) + (n_clust, -1, -1)))
        )
        msa_emb = self.linear_msa_m(msa) + tf_m

        return msa_emb, pair_emb

class SPEmbedder(nn.Module):
    """
    Embeds the spin label input in 100 bins format

    """

    def __init__(
        self,
        c_m: int,
        c_z: int,
        distograms: bool,
        **kwargs,
    ):
        """
        Args:
            c_m:
                MSA channel dimension
            c_z:
                Pair embedding channel dimension
        """
        super(SPEmbedder, self).__init__()

        self.c_m = c_m
        self.c_z = c_z

        if distograms:
            # print("@@@@@@@@@@@@@@@@", 100, self.c_z)
            self.linear = Linear(100, self.c_z)
        else:
            self.linear = Linear(1, self.c_z)
            

    def forward(
        self,
        x: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        return self.linear(x)


class RecyclingEmbedder(nn.Module):
    """
    Embeds the output of an iteration of the model for recycling.

    Implements Algorithm 32.
    """

    def __init__(
        self,
        c_m: int,
        c_z: int,
        min_bin: float,
        max_bin: float,
        no_bins: int,
        inf: float = 1e8,
        **kwargs,
    ):
        """
        Args:
            c_m:
                MSA channel dimension
            c_z:
                Pair embedding channel dimension
            min_bin:
                Smallest distogram bin (Angstroms)
            max_bin:
                Largest distogram bin (Angstroms)
            no_bins:
                Number of distogram bins
        """
        super(RecyclingEmbedder, self).__init__()

        self.c_m = c_m
        self.c_z = c_z
        self.min_bin = min_bin
        self.max_bin = max_bin
        self.no_bins = no_bins
        self.inf = inf

        self.bins = None

        self.linear = Linear(self.no_bins, self.c_z)
        self.layer_norm_m = nn.LayerNorm(self.c_m)
        self.layer_norm_z = nn.LayerNorm(self.c_z)

    def forward(
        self,
        m: torch.Tensor,
        z: torch.Tensor,
        x: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Args:
            m:
                First row of the MSA embedding. [*, N_res, C_m]
            z:
                [*, N_res, N_res, C_z] pair embedding
            x:
                [*, N_res, 3] predicted C_beta coordinates
        Returns:
            m:
                [*, N_res, C_m] MSA embedding update
            z:
                [*, N_res, N_res, C_z] pair embedding update
        """
        if self.bins is None:
            self.bins = torch.linspace(
                self.min_bin,
                self.max_bin,
                self.no_bins,
                dtype=x.dtype,
                device=x.device,
                requires_grad=False,
            )

        # [*, N, C_m]
        m_update = self.layer_norm_m(m)

        # This squared method might become problematic in FP16 mode.
        # I'm using it because my homegrown method had a stubborn discrepancy I
        # couldn't find in time.
        squared_bins = self.bins ** 2
        upper = torch.cat(
            [squared_bins[1:], squared_bins.new_tensor([self.inf])], dim=-1
        )
        d = torch.sum(
            (x[..., None, :] - x[..., None, :, :]) ** 2, dim=-1, keepdims=True
        )

        # [*, N, N, no_bins]
        d = ((d > squared_bins) * (d < upper)).type(x.dtype)

        # [*, N, N, C_z]
        d = self.linear(d)
        z_update = d + self.layer_norm_z(z)

        return m_update, z_update


class TemplateAngleEmbedder(nn.Module):
    """
    Embeds the "template_angle_feat" feature.

    Implements Algorithm 2, line 7.
    """

    def __init__(
        self,
        c_in: int,
        c_out: int,
        **kwargs,
    ):
        """
        Args:
            c_in:
                Final dimension of "template_angle_feat"
            c_out:
                Output channel dimension
        """
        super(TemplateAngleEmbedder, self).__init__()

        self.c_out = c_out
        self.c_in = c_in

        self.linear_1 = Linear(self.c_in, self.c_out, init="relu")
        self.relu = nn.ReLU()
        self.linear_2 = Linear(self.c_out, self.c_out, init="relu")

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: [*, N_templ, N_res, c_in] "template_angle_feat" features
        Returns:
            x: [*, N_templ, N_res, C_out] embedding
        """
        x = self.linear_1(x)
        x = self.relu(x)
        x = self.linear_2(x)

        return x


class TemplatePairEmbedder(nn.Module):
    """
    Embeds "template_pair_feat" features.

    Implements Algorithm 2, line 9.
    """

    def __init__(
        self,
        c_in: int,
        c_out: int,
        **kwargs,
    ):
        """
        Args:
            c_in:

            c_out:
                Output channel dimension
        """
        super(TemplatePairEmbedder, self).__init__()

        self.c_in = c_in
        self.c_out = c_out

        # Despite there being no relu nearby, the source uses that initializer
        self.linear = Linear(self.c_in, self.c_out, init="relu")

    def forward(
        self,
        x: torch.Tensor,
    ) -> torch.Tensor:
        """
        Args:
            x:
                [*, C_in] input tensor
        Returns:
            [*, C_out] output tensor
        """
        x = self.linear(x)

        return x


class ExtraMSAEmbedder(nn.Module):
    """
    Embeds unclustered MSA sequences.

    Implements Algorithm 2, line 15
    """

    def __init__(
        self,
        c_in: int,
        c_out: int,
        **kwargs,
    ):
        """
        Args:
            c_in:
                Input channel dimension
            c_out:
                Output channel dimension
        """
        super(ExtraMSAEmbedder, self).__init__()

        self.c_in = c_in
        self.c_out = c_out

        self.linear = Linear(self.c_in, self.c_out)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x:
                [*, N_extra_seq, N_res, C_in] "extra_msa_feat" features
        Returns:
            [*, N_extra_seq, N_res, C_out] embedding
        """
        x = self.linear(x)

        return x
