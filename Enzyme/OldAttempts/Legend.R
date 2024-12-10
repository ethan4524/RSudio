# ========================= Legend ========================= #
# Model Variables:
# E  - Total Extrase concentration (enzyme that converts Suproblatide to Problatide)
# S  - Concentration of Suproblatide (precursor to Problatide)
# P  - Concentration of Problatide (active ligand)
# PR - Concentration of the Problatide-GF-R receptor-ligand complex
# R  - Concentration of GF-R receptors
# Et - Active concentration of the enzyme-substrate complex (E bound to S)
# ES - Enzyme-Substrate complex (active form of Extrase interacting with Suproblatide)
# PRss - Steady state concentration of the receptor-ligand complex
# Km  - Michaelis-Menten constant for enzyme binding with Suproblatide
# Kd  - Dissociation constant for the Problatide-GF-R receptor-ligand complex
# Kcat - Catalysis rate constant for the conversion of Suproblatide to Problatide
# Kdegp - Problatide degradation/elimination rate constant
# Kdege - Extrase degradation/elimination rate constant
# Prode - Rate constant for Extrase production
# n   - Feedback parameter governing how the receptor-ligand complex affects Extrase production

# Differential Equations:
# dE   - Rate of change of Extrase concentration (total enzyme)
# dP   - Rate of change of Problatide concentration
# dPR  - Rate of change of Problatide-GF-R complex concentration

# 1. dE = (Prode * (PRss / PR)^n) - Kdege * E
#    - Describes the dynamics of Extrase, considering production of Extrase influenced by 
#      the fraction of receptor occupied by Problatide (feedback effect) and degradation of Extrase.

# 2. dP = Kcat * ES - Kdegp * P
#    - Describes the dynamics of Problatide, where Problatide is produced via catalysis by Extrase 
#      (Enzyme-Substrate complex) and degraded over time.

# 3. dPR = (P * R) / (P + Kd) - PR
#    - Describes the dynamics of the receptor-ligand complex (Problatide-GF-R), where the 
#      binding of Problatide to GF-R receptors forms the complex, and the complex decays over time.

# ========================= End Legend ======================= #
