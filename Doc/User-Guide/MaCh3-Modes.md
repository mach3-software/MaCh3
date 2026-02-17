# MaCh3 Modes {#mach3-modes}
Different neutrino generators have different naming conventions for modes in addition experiment may not want to use the same modes as the generator. This is what MaCh3 modes are used for.

This is an example of a config file to initialise modes from config.
```
#########################################
Title: "T2K Interaction modes"

GeneratorName: "NEUT"

MaCh3Modes: [
  "CCQE",
  "CC1pipm",
  "2p2h",
  ]

CCQE:
  Name: "CCQE"
  GeneratorMaping: [1]
  # 1 = CCQE

CC1pipm:
  Name: "CC 1#pi^{#pm}"
  GeneratorMaping: [11, 13]
  # 11 = CC 1pi+ 1p
  # 13 = CC 1pi+ 1n

2p2h:
  Name: "2p2h"
  GeneratorMaping: [2]
  # 2 = (Nieves) 2p2h aka MEC
```
This is what you will get:

<img width="800" src="https://github.com/mach3-software/MaCh3/assets/45295406/0af5a338-d6ef-4bad-ade0-e2bc58c3548c">


For example, NEUT has modes which are undefined, therefore all such modes are being treated `UNKNOWN_BAD`, this will not throw an error, as there may be cases when it is useful however it is worth investigating.

## How to use (Examples)
To access the actual mode you can use the following code
```
MaCh3Modes_t kMaCh3_CCQE = Modes->GetMode("CCQE");
```
To perform loop over all modes:
```
for (_int_ j = 0; j < Modes->GetNModes()+1; j++)
{
    Modes->GetMaCh3ModeName(j);
}
```

## Unknown Category
If for some reason you have gaps in your mode definition or will try to access more modes than have been defined you will get the uknonw category. This is to prevent segfaults, also there are cases when your potential number of modes is too large and having a category where all garbage categories are assigned is useful.

