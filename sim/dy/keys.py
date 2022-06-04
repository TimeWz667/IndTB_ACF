# Dim 0: TB
U = 0

# Active TB
Asym_Sn_DS = 1
Asym_Sp_DS = 2
Asym_Sn_DR = 3
Asym_Sp_DR = 4

Sym_Sn_DS = 5
Sym_Sp_DS = 6
Sym_Sn_DR = 7
Sym_Sp_DR = 8

ExSym_Sn_DS = 9
ExSym_Sp_DS = 10
ExSym_Sn_DR = 11
ExSym_Sp_DR = 12

# Tx
Txf_Pub_Sn_DS = 13
Txf_Pub_Sp_DS = 14
Txf_Pub_Sn_DR = 15
Txf_Pub_Sp_DR = 16

Txf_Pri_Sn_DS = 17
Txf_Pri_Sp_DS = 18
Txf_Pri_Sn_DR = 19
Txf_Pri_Sp_DR = 20

Txs_Pub_Sn_DS = 21
Txs_Pub_Sp_DS = 22
Txs_Pub_Sn_DR = 23
Txs_Pub_Sp_DR = 24

Txs_Pri_Sn_DS = 25
Txs_Pri_Sp_DS = 26
Txs_Pri_Sn_DR = 27
Txs_Pri_Sp_DR = 28

FLat_DS = 29
FLat_DR = 30
SLat_DS = 31
SLat_DR = 32
RLow_DS = 33
RLow_DR = 34
RHigh_DS = 35
RHigh_DR = 36
RSt_DS = 37
RSt_DR = 38

N_State_TB = 39


# Meta group
Asym = [Asym_Sn_DS, Asym_Sp_DS, Asym_Sn_DR, Asym_Sp_DR]
Sym = [Sym_Sn_DS, Sym_Sp_DS, Sym_Sn_DR, Sym_Sp_DR]
ExSym = [ExSym_Sn_DS, ExSym_Sp_DS, ExSym_Sn_DR, ExSym_Sp_DR]

Asym_Sn = [Asym_Sn_DS, Asym_Sn_DR]
Sym_Sn = [Sym_Sn_DS, Sym_Sn_DR]
ExSym_Sn = [ExSym_Sn_DS, ExSym_Sn_DR]

Asym_Sp = [Asym_Sp_DS, Asym_Sp_DR]
Sym_Sp = [Sym_Sp_DS, Sym_Sp_DR]
ExSym_Sp = [ExSym_Sp_DS, ExSym_Sp_DR]
# Pub = [Pub_Sn_DS, Pub_Sp_DS, Pub_Sn_DR, Pub_Sp_DR]
# Pri = [Pri_Sn_DS, Pri_Sp_DS, Pri_Sn_DR, Pri_Sp_DR]
# ACFS_Asym = [ACFS_Asym_Sn_DS, ACFS_Asym_Sp_DS, ACFS_Asym_Sn_DR, ACFS_Asym_Sp_DR]
# ACFS_Sym = [ACFS_Sym_Sn_DS, ACFS_Sym_Sp_DS, ACFS_Sym_Sn_DR, ACFS_Sym_Sp_DR]
# ACFC_Sym = [ACFC_Sym_Sn_DS, ACFC_Sym_Sp_DS, ACFC_Sym_Sn_DR, ACFC_Sym_Sp_DR]
Txf_Pub = [Txf_Pub_Sn_DS, Txf_Pub_Sp_DS, Txf_Pub_Sn_DR, Txf_Pub_Sp_DR]
Txf_Pri = [Txf_Pri_Sn_DS, Txf_Pri_Sp_DS, Txf_Pri_Sn_DR, Txf_Pri_Sp_DR]
Txs_Pub = [Txs_Pub_Sn_DS, Txs_Pub_Sp_DS, Txs_Pub_Sn_DR, Txs_Pub_Sp_DR]
Txs_Pri = [Txs_Pri_Sn_DS, Txs_Pri_Sp_DS, Txs_Pri_Sn_DR, Txs_Pri_Sp_DR]

FLat = [FLat_DS, FLat_DR]
SLat = [SLat_DS, SLat_DR]
RLow = [RLow_DS, RLow_DR]
RHigh = [RHigh_DS, RHigh_DR]
RSt = [RSt_DS, RSt_DR]


# Infectious
Infectious_Sn_DS = [x[0] for x in [Asym, Sym, ExSym]]
Infectious_Sp_DS = [x[1] for x in [Asym, Sym, ExSym]]
Infectious_Sn_DR = [x[2] for x in [Asym, Sym, ExSym]]
Infectious_Sp_DR = [x[3] for x in [Asym, Sym, ExSym]]
# Infectious_Sn_DS = [x[0] for x in [Asym, Sym, ExSym, Pub, Pri, ACFS_Asym, ACFS_Sym, ACFC_Sym]]
# Infectious_Sp_DS = [x[1] for x in [Asym, Sym, ExSym, Pub, Pri, ACFS_Asym, ACFS_Sym, ACFC_Sym]]
# Infectious_Sn_DR = [x[2] for x in [Asym, Sym, ExSym, Pub, Pri, ACFS_Asym, ACFS_Sym, ACFC_Sym]]
# Infectious_Sp_DR = [x[3] for x in [Asym, Sym, ExSym, Pub, Pri, ACFS_Asym, ACFS_Sym, ACFC_Sym]]

Infectious = Infectious_Sn_DS + Infectious_Sp_DS + Infectious_Sn_DR + Infectious_Sp_DR

Infectious_Sn_DS.sort()
Infectious_Sp_DS.sort()
Infectious_Sn_DR.sort()
Infectious_Sp_DR.sort()
Infectious.sort()


# Healthcare
# PCF = Pub + Pri
# PCF.sort()
#
# ACF = ACFS_Asym + ACFS_Sym + ACFC_Sym
# ACF.sort()


# LTBI
LTBI = FLat + SLat + RLow + RHigh + RSt
LTBI.sort()

Sub_Resist = Sub_DS, Sub_DR = [0, 1]  # q
Sub_Smear = Sub_Sn, Sub_Sp = [0, 1]  # r

Dim_qr = (len(Sub_Resist), len(Sub_Smear))


# Second dim
N_State_Strata = 2
RiskLo, RiskHi = 0, 1
Tag_Strata = ['RiskLo', 'RiskHi']
