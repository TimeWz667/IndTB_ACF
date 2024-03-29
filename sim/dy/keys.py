# Dim 0: TB
U = 0

# Active TB
Asym_DS = 1
Asym_DR = 2

Sym_DS = 3
Sym_DR = 4

ExSym_DS = 5
ExSym_DR = 6

# Tx
Txf_Pub_DS = 7
Txf_Pub_DR = 8

Txf_Pri_DS = 9
Txf_Pri_DR = 10

Txs_Pub_DS = 11
Txs_Pub_DR = 12

Txs_Pri_DS = 13
Txs_Pri_DR = 14

FLat_DS = 15
FLat_DR = 16
SLat_DS = 17
SLat_DR = 18
RLow_DS = 19
RLow_DR = 20
RHigh_DS = 21
RHigh_DR = 22
RSt_DS = 23
RSt_DR = 24

UTPT = 25

FLatTPT_DS = 26
FLatTPT_DR = 27
SLatTPT_DS = 28
SLatTPT_DR = 29
RLowTPT_DS = 30
RLowTPT_DR = 31
RHighTPT_DS = 32
RHighTPT_DR = 33
RStTPT_DS = 34
RStTPT_DR = 35


N_State_TB = 36


# Meta group
Asym = [Asym_DS, Asym_DR]
Sym = [Sym_DS, Sym_DR]
ExSym = [ExSym_DS, ExSym_DR]

Txf_Pub = [Txf_Pub_DS, Txf_Pub_DR]
Txf_Pri = [Txf_Pri_DS, Txf_Pri_DR]
Txs_Pub = [Txs_Pub_DS, Txs_Pub_DR]
Txs_Pri = [Txs_Pri_DS, Txs_Pri_DR]
Tx_DS = [Txf_Pub_DS, Txf_Pri_DS, Txs_Pub_DS, Txs_Pri_DS]
Tx_DR = [Txf_Pub_DR, Txf_Pri_DR, Txs_Pub_DR, Txs_Pri_DR]

FLat = [FLat_DS, FLat_DR]
SLat = [SLat_DS, SLat_DR]
RLow = [RLow_DS, RLow_DR]
RHigh = [RHigh_DS, RHigh_DR]
RSt = [RSt_DS, RSt_DR]

FLatTPT = [FLatTPT_DS, FLatTPT_DR]
SLatTPT = [SLatTPT_DS, SLatTPT_DR]
RLowTPT = [RLowTPT_DS, RLowTPT_DR]
RHighTPT = [RHighTPT_DS, RHighTPT_DR]
RStTPT = [RStTPT_DS, RStTPT_DR]


# Infectious
Infectious_DS = [x[0] for x in [Asym, Sym, ExSym]]
Infectious_DR = [x[1] for x in [Asym, Sym, ExSym]]

Infectious = Infectious_DS + Infectious_DR

Infectious_DS.sort()
Infectious_DR.sort()
Infectious.sort()


# LTBI
LTBI0 = FLat + SLat + RLow + RHigh + RSt
LTBI0.sort()

LTBI_TPT = FLatTPT + SLatTPT + RLowTPT + RHighTPT + RStTPT
LTBI_TPT.sort()

LTBI = LTBI0 + LTBI_TPT
LTBI.sort()


# Second dim
N_State_Strata = 2
RiskLo, RiskHi = 0, 1
Tag_Strata = ['RiskLo', 'RiskHi']

N_Aux = 2
IncN = 0
MorN = 1
