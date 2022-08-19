
#include "yaml-cpp/yaml.h"

#include <iostream>

int main(int argc, char const *argv[]) {

  std::string doc_str = R"(
keys: [name generated prior lowerbound upperbound error renorm type mode elements nupdg q2_range etru_range splineind detid sk_spline_name nd_spline_name stepscale sk_mode flatprior]
parameters:
  -
    name: MAQE
    generated: 1.21
    prior: 1.03
    lowerbound: 0
    upperbound: 9999
    error: 0.06
    renorm: 0
    type: spline
    splineind: 0
    detid: 985
    sk_spline_name: MaCCQE
    nd_spline_name: MAQEGraph
    stepscale: 7.5
    sk_mode: 0
  -
    name: Q2_norm_5
    generated: 1.0
    prior: 1.0
    lowerbound: 0
    upperbound: 9999
    error: 0.11
    renorm: 0
    type: norm
    detid: 985
    stepscale: 3.0
    mode: 0
    elements: -1 12 16
    q2_range: 0.25 0.5
  -
    name: Q2_norm_6
    generated: 1.0
    prior: 1.0
    lowerbound: 0
    upperbound: 9999
    error: 0.18
    renorm: 0
    type: norm
    detid: 985
    stepscale: 3.0
    mode: 0
    elements: -1 12 16
    q2_range: 0.5 1.0
  -
    name: Q2_norm_7
    generated: 1.0
    prior: 1.0
    lowerbound: 0
    upperbound: 9999
    error: 0.40
    renorm: 0
    type: norm
    detid: 985
    stepscale: 3.0
    mode: 0
    elements: -1 12 16
    q2_range: 1.0 -999
  -
    name: PShell_MF_Norm_C
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 9999
    error: 0.2
    renorm: 0
    type: spline
    splineind: 1
    detid: 1
    nd_spline_name: PShell_MF_Norm_CGraph
    stepscale: 4.2
    sk_mode: 0
  -
    name: SShell_MF_Norm_C
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 9999
    error: 0.45
    renorm: 0
    type: spline
    splineind: 2
    detid: 1
    nd_spline_name: SShell_MF_Norm_CGraph
    stepscale: 4.2
    sk_mode: 0
  -
    name: SRC_Norm_C
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 2
    renorm: 0
    type: spline
    splineind: 3
    detid: 1
    sk_spline_name: SFSRCNormC
    nd_spline_name: SRC_Norm_CGraph
    stepscale: 1.2
    sk_mode: 0
  -
    name: PShell_MF_PMissShape_C
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1
    renorm: 0
    type: spline
    splineind: 4
    detid: 1
    nd_spline_name: PShell_MF_PMissShape_CGraph
    stepscale: 2.0
    sk_mode: 0
  -
    name: SShell_MF_PMissShape_C
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1
    renorm: 0
    type: spline
    splineind: 5
    detid: 1
    nd_spline_name: SShell_MF_PMissShape_CGraph
    stepscale: 3.0
    sk_mode: 0
  -
    name: P1_2Shell_MF_Norm_O
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 9999
    error: 0.2
    renorm: 0
    type: spline
    splineind: 6
    detid: 985
    sk_spline_name: SFP12ShellMeanFNormO
    nd_spline_name: P1_2Shell_MF_Norm_OGraph
    stepscale: 3.0
    sk_mode: 0
  -
    name: P3_2Shell_MF_Norm_O
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 9999
    error: 0.45
    renorm: 0
    type: spline
    splineind: 7
    detid: 985
    sk_spline_name: SFP32ShellMeanFNormO
    nd_spline_name: P3_2Shell_MF_Norm_OGraph
    stepscale: 3.0
    sk_mode: 0
  -
    name: SShell_MF_Norm_O
    sk_spline_name: SFSShellMeanFNormO
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 9999
    error: 0.75
    renorm: 0
    type: spline
    splineind: 8
    detid: 985
    nd_spline_name: SShell_MF_Norm_OGraph
    stepscale: 3.0
    sk_mode: 0
  -
    name: SRC_Norm_O
    sk_spline_name: SFSRCNormO
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 2
    renorm: 0
    type: spline
    splineind: 9
    detid: 985
    nd_spline_name: SRC_Norm_OGraph
    stepscale: 1.0
    sk_mode: 0
  -
    name: P1_2Shell_MF_PMissShape_O
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1
    renorm: 0
    type: spline
    splineind: 10
    detid: 985
    sk_spline_name: SFP12ShellMeanFPMissShapeO
    nd_spline_name: P1_2Shell_MF_PMissShape_OGraph
    stepscale: 1.0
    sk_mode: 0
  -
    name: P3_2Shell_MF_PMissShape_O
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1
    renorm: 0
    type: spline
    splineind: 11
    detid: 985
    sk_spline_name: SFP32ShellMeanFPMissShapeO
    nd_spline_name: P3_2Shell_MF_PMissShape_OGraph
    stepscale: 2.0
    sk_mode: 0
  -
    name: SShell_MF_PMissShape_O
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1
    renorm: 0
    type: spline
    splineind: 12
    detid: 985
    sk_spline_name: SFSShellMeanFPMissShapeO
    nd_spline_name: SShell_MF_PMissShape_OGraph
    stepscale: 2.0
    sk_mode: 0
  -
    name: Pauli_Blocking_C_nu
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1
    renorm: 0
    type: spline
    splineind: 13
    detid: 1
    nd_spline_name: Pauli_Blocking_C12_nuGraph
    stepscale: 3.5
    correlation:
      Pauli_Blocking_C_nubar: 0.8
    sk_mode: 0
  -
    name: Pauli_Blocking_O_nu
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1
    renorm: 0
    type: spline
    splineind: 14
    detid: 985
    sk_spline_name: SFPBTwkDialHybridO16nu
    nd_spline_name: Pauli_Blocking_O16_nuGraph
    stepscale: 3.5
    correlation:
      Pauli_Blocking_O_nubar: 0.8
    sk_mode: 0
  -
    name: Pauli_Blocking_C_nubar
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1
    renorm: 0
    type: spline
    splineind: 15
    detid: 1
    nd_spline_name: Pauli_Blocking_C12_nubarGraph
    stepscale: 3.5
    correlation:
      Pauli_Blocking_C_nu: 0.8
    sk_mode: 0
  -
    name: Pauli_Blocking_O_nubar
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1
    renorm: 0
    type: spline
    splineind: 16
    detid: 985
    sk_spline_name: SFPBTwkDialHybridO16nubar
    nd_spline_name: Pauli_Blocking_O16_nubarGraph
    stepscale: 3.5
    correlation:
      Pauli_Blocking_O_nu: 0.8
    sk_mode: 0
  -
    name: Optical_Potential_C
    generated: 0
    prior: 0
    lowerbound: 0
    upperbound: 1
    error: 1
    renorm: 0
    type: spline
    splineind: 17
    detid: 1
    nd_spline_name: Optical_Potential_C12Graph
    stepscale: 1.0
    sk_mode: 0
    flatprior: 1
  -
    name: Optical_Potential_O
    generated: 0
    prior: 0
    lowerbound: 0
    upperbound: 1
    error: 1
    renorm: 0
    type: spline
    splineind: 18
    detid: 985
    sk_spline_name: SFOptPotTwkDialO16
    nd_spline_name: Optical_Potential_O16Graph
    stepscale: 1.0
    sk_mode: 0
    flatprior: 1
  -
    name: 2p2h_norm_nu
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 1
    renorm: 0
    type: norm
    detid: 985
    stepscale: 2.0
    mode: 9
    element: -1 12 16
    nupdg: 12 14 16
    flatprior: 1
  -
    name: 2p2h_norm_nubar
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 1
    renorm: 0
    type: norm
    detid: 985
    stepscale: 2.0
    mode: 9
    element: -1 12 16
    nupdg: -12 -14 -16
    flatprior: 1
  -
    name: 2p2h_normCtoO
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 0.2
    renorm: 0
    type: norm
    detid: 985
    stepscale: 4.0
    mode: 9
    element: 16
  -
    name: 2p2h_Edep_lowEnu
    generated: 1.0
    prior: 1.00
    lowerbound: 0
    upperbound: 1.0
    error: 1
    renorm: 0
    type: spline
    splineind: 19
    detid: 985
    sk_spline_name: 2p2hedeplowenu
    nd_spline_name: MEC_lowEnuGraph
    stepscale: 1.0
    sk_mode: 9
    flatprior: 1
  -
    name: 2p2h_Edep_highEnu
    generated: 1.0
    prior: 1.0
    lowerbound: 0
    upperbound: 1.0
    error: 1
    renorm: 0
    type: spline
    splineind: 20
    detid: 985
    sk_spline_name: 2p2hedephienu
    nd_spline_name: MEC_highEnuGraph
    stepscale: 1.0
    sk_mode: 9
    flatprior: 1
  -
    name: 2p2h_Edep_lowEnubar
    generated: 1.0
    prior: 1.0
    lowerbound: 0
    upperbound: 1.0
    error: 1
    renorm: 0
    type: spline
    splineind: 21
    detid: 985
    sk_spline_name: 2p2hedeplowenubar
    nd_spline_name: MEC_lowEnubarGraph
    stepscale: 1.0
    sk_mode: 9
    flatprior: 1
  -
    name: 2p2h_Edep_highEnubar
    generated: 1.0
    prior: 1.0
    lowerbound: 0
    upperbound: 1.0
    error: 1
    renorm: 0
    type: spline
    splineind: 22
    detid: 985
    sk_spline_name: 2p2hedephienubar
    nd_spline_name: MEC_highEnubarGraph
    stepscale: 1.0
    sk_mode: 9
    flatprior: 1
  -
    name: PNNN_Shape
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 1
    error: 0.33
    renorm: 0
    type: spline
    splineind: 23
    detid: 985
    sk_spline_name: MECTwkDialPNNNShape
    nd_spline_name: PNNN_ShapeGraph
    stepscale: 2.2
    sk_mode: 9
  -
    name: 2p2h_shape_C_np
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 1
    error: 3.0
    renorm: 0
    type: spline
    splineind: 24
    detid: 1
    nd_spline_name: PDD_C_npGraph
    stepscale: 0.6
    correlation:
      2p2h_shape_O_np: 0.3
    sk_mode: 9
  -
    name: 2p2h_shape_C_NN
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 1
    error: 3.0
    renorm: 0
    type: spline
    splineind: 25
    detid: 1
    nd_spline_name: PDD_C_NNGraph
    stepscale: 0.6
    correlation:
      2p2h_shape_O_NN: 0.3
    sk_mode: 9
  -
    name: 2p2h_shape_O_np
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 1
    error: 3.0
    renorm: 0
    type: spline
    splineind: 26
    detid: 985
    sk_spline_name: MECTwkDialPDDWeightO16np
    nd_spline_name: PDD_O_npGraph
    stepscale: 0.6
    correlation:
      2p2h_shape_C_np: 0.3
    sk_mode: 9
  -
    name: 2p2h_shape_O_NN
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 1
    error: 3.0
    renorm: 0
    type: spline
    splineind: 27
    detid: 985
    sk_spline_name: MECTwkDialPDDWeightO16NN
    nd_spline_name: PDD_O_NNGraph
    stepscale: 0.6
    correlation:
      2p2h_shape_C_NN: 0.3
    sk_mode: 9
  -
    name: CA5
    generated: 1.01
    prior: 1.06
    lowerbound: 0
    upperbound: 9999
    error: 0.1
    renorm: 0
    type: spline
    splineind: 28
    detid: 985
    sk_spline_name: CA5RES
    nd_spline_name: CA5Graph
    stepscale: 3.0
    correlation:
      MARES: -0.11
      ISO_BKG: -0.03
    sk_mode: 1 5 6 14
  -
    name: MARES
    generated: 0.95
    prior: 0.91
    lowerbound: 0
    upperbound: 9999
    error: 0.1
    renorm: 0
    type: spline
    splineind: 29
    detid: 985
    sk_spline_name: MaRES
    nd_spline_name: MARESGraph
    stepscale: 4.0
    correlation:
      CA5: -0.11
    sk_mode: 1 5 6 14
  -
    name: ISO_BKG_LowPPi
    generated: 1.3
    prior: 1.3
    lowerbound: 0
    upperbound: 9999
    error: 1.30
    renorm: 0
    type: spline
    splineind: 30
    detid: 984
    sk_spline_name: bgsclbar
    nd_spline_name: BgSclRes_lowPPiGraph
    stepscale: 4.0
    sk_mode: 1 5 6 14
  -
    name: ISO_BKG
    generated: 1.3
    prior: 1.21
    lowerbound: 0
    upperbound: 9999
    error: 0.27
    renorm: 0
    type: spline
    splineind: 31
    detid: 985
    sk_spline_name: bgscl
    nd_spline_name: BgSclResGraph
    stepscale: 4.0
    correlation:
      CA5: -0.03
    sk_mode: 1 5 6 14
  -
    name: RES_Eb_C_numu
    generated: 0
    prior: 25
    lowerbound: 0
    upperbound: 50
    error: 25.0
    renorm: 0
    type: spline
    splineind: 32
    detid: 1
    nd_spline_name: RES_Eb_C_numuGraph
    stepscale: 1.5
    sk_mode: 1 5 6 14
  -
    name: RES_Eb_O_numu
    generated: 0
    prior: 25
    lowerbound: 0
    upperbound: 50
    error: 25.0
    renorm: 0
    type: spline
    splineind: 33
    detid: 985
    sk_spline_name: RESEbOnumu
    nd_spline_name: RES_Eb_O_numuGraph
    stepscale: 1.5
    sk_mode: 1 5 6 14
  -
    name: RES_Eb_C_numubar
    generated: 0
    prior: 25
    lowerbound: 0
    upperbound: 50
    error: 25.0
    renorm: 0
    type: spline
    splineind: 34
    detid: 1
    nd_spline_name: RES_Eb_C_numubarGraph
    stepscale: 1.5
    sk_mode: 1 5 6 14
  -
    name: RES_Eb_O_numubar
    generated: 0
    prior: 25
    lowerbound: 0
    upperbound: 50
    error: 25.0
    renorm: 0
    type: spline
    splineind: 35
    detid: 985
    sk_spline_name: RESEbOnumubar
    nd_spline_name: RES_Eb_O_numubarGraph
    stepscale: 1.5
    sk_mode: 1 5 6 14
  -
    name: RS_Delta_Decay
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 1.0
    error: 1.0
    renorm: 0
    type: spline
    splineind: 36
    detid: 985
    sk_spline_name: MDLSPiEj
    nd_spline_name: R_S_Delta_DecayGraph
    stepscale: 1.2
    sk_mode: 1 5 6 14
    flatprior: 1
  -
    name: SPP_Pi0Norm_numu
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 0.3
    renorm: 0
    type: norm
    detid: 985
    stepscale: 4.0
    mode: 14
    nupdg: 14
  -
    name: SPP_Pi0Norm_numubar
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 0.3
    renorm: 0
    type: norm
    detid: 985
    stepscale: 4.0
    mode: 14
    nupdg: -14
  -
    name: FEFQE
    generated: 1.069
    prior: 1.069
    lowerbound: 0
    upperbound: 9999
    error: 0.313
    renorm: 0
    type: spline
    splineind: 37
    detid: 985
    sk_spline_name: PionFSIQELow
    nd_spline_name: PionFSI_QELowMomProbGraph
    stepscale: 2.0
    correlation:
      FEFABS: -0.27
      FEFCX: -0.20
      FEFINEL: 0.03
      FEFQEH: -0.03
    sk_mode: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  -
    name: FEFQEH
    generated: 1.824
    prior: 1.824
    lowerbound: 0
    upperbound: 9999
    error: 0.859
    renorm: 0
    type: spline
    splineind: 38
    detid: 985
    sk_spline_name: PionFSIQEHigh
    nd_spline_name: PionFSI_QEHighMomProbGraph
    stepscale: 2.0
    correlation:
      FEFQE: -0.03
      FEFABS: 0.03
      FEFCX: 0.01
      FEFINEL: -0.53
    sk_mode: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  -
    name: FEFINEL
    generated: 1.002
    prior: 1.002
    lowerbound: 0
    upperbound: 9999
    error: 1.101
    renorm: 0
    type: spline
    splineind: 39
    detid: 985
    sk_spline_name: PionFSIInel
    nd_spline_name: PionFSI_InelProbGraph
    stepscale: 2.0
    correlation:
      FEFQE: 0.03
      FEFABS: -0.03
      FEFCX: -0.01
      FEFQEH: -0.53
    sk_mode: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  -
    name: FEFABS
    generated: 1.404
    prior: 1.404
    lowerbound: 0
    upperbound: 9999
    error: 0.432
    renorm: 0
    type: spline
    splineind: 40
    detid: 985
    sk_spline_name: PionFSIAbs
    nd_spline_name: PionFSI_AbsProbGraph
    stepscale: 2.0
    correlation:
      FEFQE: -0.27
      FEFCX: -0.31
      FEFINEL: -0.03
      FEFQEH: 0.03
    sk_mode: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  -
    name: FEFCX
    generated: 0.697
    prior: 0.697
    lowerbound: 0
    upperbound: 9999
    error: 0.305
    renorm: 0
    type: spline
    splineind: 41
    detid: 985
    sk_spline_name: PionFSICExLow
    nd_spline_name: PionFSI_CExLowMomProbGraph
    stepscale: 3.0
    correlation:
      FEFQE: -0.20
      FEFABS: -0.31
      FEFINEL: -0.01
      FEFQEH: 0.01
    sk_mode: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  -
    name: FEFCXH
    generated: 1.8
    prior: 1.8
    lowerbound: 0
    upperbound: 9999
    error: 0.288
    renorm: 0
    type: spline
    splineind: 42
    detid: 985
    sk_spline_name: PionFSICExHigh
    nd_spline_name: PionFSI_CExHighMomProbGraph
    stepscale: 3.2
    sk_mode: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
  -
    name: Nucleon_FSI
    generated: 0
    prior: 0
    lowerbound: -1
    upperbound: 1
    error: 0.3
    renorm: 0
    type: spline
    splineind: 43
    detid: 985
    sk_spline_name: TwkDialFateNucleonFSI
    nd_spline_name: Nucleon_Fate_FSIGraph
    stepscale: 3.0
    sk_mode: 0 1 2 9 14
  -
    name: CC_Coh_C
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 0.3
    renorm: 0
    type: norm
    detid: 1
    stepscale: 4.0
    correlation:
      CC_Coh_O: 0.9999
    mode: 2
    element: 12
  -
    name: CC_Coh_O
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 0.3
    renorm: 0
    type: norm
    detid: 985
    stepscale: 4.0
    correlation:
      CC_Coh_C: 0.9999
    mode: 2
    element: 16
  -
    name: MPi_Multi_TotXSec
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1.0
    renorm: 0
    type: spline
    splineind: 44
    detid: 985
    sk_spline_name: MPiTotXsec
    nd_spline_name: MultiPi_Multiplicity_TotXSecGraph
    stepscale: 3.0
    sk_mode: 3
  -
    name: MPi_BY_Vector
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1.0
    renorm: 0
    type: spline
    splineind: 45
    detid: 985
    sk_spline_name: MultiPiBYVector
    nd_spline_name: MultiPi_BY_VectorGraph
    stepscale: 3.5
    sk_mode: 3
  -
    name: MPi_BY_Axial
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1.0
    renorm: 0
    type: spline
    splineind: 46
    detid: 985
    sk_spline_name: MultiPiBYAxial
    nd_spline_name: MultiPi_BY_AxialGraph
    stepscale: 3.0
    sk_mode: 3
  -
    name: MPi_Multi_Shape
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1.0
    renorm: 0
    type: spline
    splineind: 47
    detid: 985
    sk_spline_name: MultiPiMultiplicityShape
    nd_spline_name: MultiPi_Multiplicity_ShapeGraph
    stepscale: 3.0
    sk_mode: 3
  -
    name: CC_BY_DIS
    generated: 0
    prior: 0
    lowerbound: -9999
    upperbound: 9999
    error: 1.0
    renorm: 0
    type: spline
    splineind: 48
    detid: 985
    sk_spline_name: DISBY
    nd_spline_name: DIS_BY_corrGraph
    stepscale: 3.0
    sk_mode: 4
  -
    name: CC_DIS_MultPi_Norm_Nu
    generated: 1.0
    prior: 1.0
    lowerbound: 0
    upperbound: 9999
    error: 0.035
    renorm: 0
    type: norm
    detid: 985
    stepscale: 5.0
    mode: 3 4
    element: 12 16
    nupdg: 12 14 16
    correlation:
      CC_DIS_MultPi_Norm_Nu: 0.999999
  -
    name: CC_DIS_MultPi_Norm_Nubar
    generated: 1.0
    prior: 1.0
    lowerbound: 0
    upperbound: 9999
    error: 0.065
    renorm: 0
    type: norm
    detid: 985
    stepscale: 5.0
    mode: 3 4
    element: 12 16
    nupdg: -12 -14 -16
    correlation:
      CC_DIS_MultPi_Norm_Nubar: 0.999999
  -
    name: CC_Misc
    generated: 1.0
    prior: 1.0
    lowerbound: 0
    upperbound: 9999
    error: 1.0
    renorm: 0
    type: norm
    detid: 985
    stepscale: 5.0
    mode: 11
    element: 12 16
  -
    name: NC_Coh
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 0.3
    renorm: 0
    type: norm
    detid: 985
    stepscale: 4.0
    mode: 7
    element:
  -
    name: NC_1gamma
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 1
    renorm: 0
    type: norm
    detid: 985
    stepscale: 4.0
    mode: 10
    element:
  -
    name: NC_other_near
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 0.3
    renorm: 0
    type: norm
    detid: 1
    stepscale: 4.0
    mode: 8 12 13
    element:
  -
    name: NC_other_far
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 0.3
    renorm: 0
    type: norm
    detid: 984
    stepscale: 4.0
    mode: 8 12 13
    element:
  -
    name: CC_norm_nu
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 0.02
    renorm: 0
    type: norm
    detid: 985
    stepscale: 4.0
    mode: 0 1 2 3 4 9 11 14
    element: 12 16 
    nupdg: 12 14 16
    etru_range: 0.3 0.6
    correlation:
      CC_norm_nubar: -0.999999
  -
    name: CC_norm_nubar
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 0.01
    renorm: 0
    type: norm
    detid: 985
    stepscale: 4.0
    mode: 0 1 2 3 4 9 11 14
    element: 12 16
    nupdg: -12 -14 -16
    etru_range: 0.3 0.6
    correlation:
      CC_norm_nu: -0.999999
  -
    name: nue_numu
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 2.82842712474619014e-02
    renorm: 0
    type: norm
    detid: 985
    stepscale: 4.0
    correlation:
      nuebar_numubar: -0.5
    mode: 0 1 2 3 4 9 11 14
    element:
    nupdg: 12
  -
    name: nuebar_numubar
    generated: 1
    prior: 1
    lowerbound: 0
    upperbound: 9999
    error: 2.82842712474619014e-02
    renorm: 0
    type: norm
    detid: 985
    stepscale: 4.0
    correlation:
      nue_numu: -0.5
    mode: 0 1 2 3 4 9 11 14
    element:
    nupdg: -12
  -
    name: EB_dial_C_nu
    generated: 0
    prior: 2
    lowerbound: -10
    upperbound: 15
    error: 6.0
    renorm: 0
    type: function
    detid: 1
    stepscale: 6.0
    correlation:
      EB_dial_C_nubar: 0.778
      EB_dial_O_nu: 0.870
      EB_dial_O_nubar: 0.653
  -
    name: EB_dial_C_nubar
    generated: 0
    prior: 0
    lowerbound: -10
    upperbound: 15
    error: 6.0
    renorm: 0
    type: function
    detid: 1
    stepscale: 6.0
    correlation:
      EB_dial_C_nu: 0.778
      EB_dial_O_nu: 0.653
      EB_dial_O_nubar: 0.870
  -
    name: EB_dial_O_nu
    generated: 0
    prior: 4
    lowerbound: -10
    upperbound: 15
    error: 6.0
    renorm: 0
    type: function
    detid: 985
    stepscale: 6.0
    correlation:
      EB_dial_C_nu: 0.870
      EB_dial_C_nubar: 0.653
      EB_dial_O_nubar: 0.778
  -
    name: EB_dial_O_nubar
    generated: 0
    prior: 0
    lowerbound: -10
    upperbound: 15
    error: 6.0
    renorm: 0
    type: function
    detid: 985
    stepscale: 6.0
    correlation:
      EB_dial_C_nu: 0.653
      EB_dial_C_nubar: 0.870
      EB_dial_O_nu: 0.778
  )";

  auto doc = YAML::Load(doc_str);

  auto keys = doc["keys"].as<std::vector<std::string>>();

  std::cout << "Keys: " << std::endl;
  for (auto k : keys) {
    std::cout << "\t" << k << std::endl;
  }
  std::cout << "\n\n" << std::endl;

  std::cout << "Parameters: " << std::endl;
  for (auto const &param : doc["parameters"]) {
    std::cout << "\tname: " << param["name"].as<std::string>() << std::endl;
    if (param["correlation"]) {
      std::cout << "\t\thas " << param["correlation"].size() << " correlations:" << std::endl;
      for (auto const &corr : param["correlation"]) {
        std::cout << "\t\t\t" << corr.first.as<std::string>() << ": " << corr.second.as<double>() << std::endl;
      }
    }
  }
}