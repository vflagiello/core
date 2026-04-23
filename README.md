# core

##############################################################################
# PREPARE VARIABLES FOR NIT CALCULATION
##############################################################################
master_nits <- master_nits %>%
mutate(
AST_U_L = ifelse(AST_U_L == 0, 0.01, AST_U_L),
ALT_U_L = ifelse(ALT_U_L == 0, 0.01, ALT_U_L),
GGT_U_L = ifelse(GGT_U_L == 0, 0.01, GGT_U_L),
bmi_capped_40 = ifelse(bmi <= 40, bmi, 40),
albumin_g_dL = albumin_g_L / 10,
total_protein_g_dL = total_protein_g_L / 10,
globulin_g_dL = total_protein_g_dL - albumin_g_dL,
ast_alt_ratio = AST_U_L / ALT_U_L,
HbA1c_percent = HbA1c_mmol_mol * 0.0915 + 2.15,
HDL_mg_dL = HDL_mmol_L * 39,
AST_ukat_L = AST_U_L / 60,
ALT_ukat_L = ALT_U_L / 60,
GGT_ukat_L = GGT_U_L / 60,
ln_ast_core = log(AST_ukat_L),
ln_alt_core = log(ALT_ukat_L),
ln_ggt_core = log(GGT_ukat_L))

##############################################################################
# CALCULATE NITs
##############################################################################
master_nits <- master_nits %>%
  mutate(
    maf5 = -11.3674 + 0.0282 * waist_c - 0.1761 * bmi + 0.0019 * (waist_c * bmi) + 2.0762 * dm2_baseline + 2.9207 * log(AST_U_L) - 0.0059 * platelets_10e9_L,
    safe = 2.97 * age + 5.99 * bmi_capped_40 + 62.85 * dm2_baseline + 154.85 * log(AST_U_L) - 58.23 * log(ALT_U_L) + 195.48 * log(globulin_g_dL) - 141.61 * log(platelets_10e9_L) - 75,
    nfs = -1.675 + 0.037 * age + 0.094 * bmi + 1.13 * dm2_baseline + 0.99 * ast_alt_ratio - 0.013 * platelets_10e9_L - 0.66 * albumin_g_dL,
    fni = exp(-10.33 + 2.54 * log(AST_U_L) + 3.86 * log(HbA1c_percent) - 1.66 * log(HDL_mg_dL)) / (1 + exp(-10.33 + 2.54 * log(AST_U_L) + 3.86 * log(HbA1c_percent) - 1.66 * log(HDL_mg_dL))),
    fib4 = (age * AST_U_L) / (platelets_10e9_L * sqrt(ALT_U_L)),
    bmi_smaf5 = ifelse(bmi < 25, 25, bmi),
    male_smaf5 = ifelse(sex == "male", 1, 0),
    smaf5 = -15.773 + 0.247 * bmi_smaf5 + 2.339 * dm2_baseline + 2.950 * log(AST_U_L) - 0.00686 * platelets_10e9_L + 0.503 * male_smaf5
  )

tp3 <- function(x, k) pmax(x - k, 0)^3
male_core <- ifelse(master_nits$sex == "male", 1, 0)

lp_core <- with(master_nits,
                -13.0001 +
                  0.1528 * age -
                  0.0016 * tp3(age, 30.65) +
                  0.0052 * tp3(age, 39.35) -
                  0.0072 * tp3(age, 46.65) +
                  0.0036 * tp3(age, 56.18) -
                  0.6865 * ln_ggt_core -
                  0.0056 * tp3(ln_ggt_core, -1.966113) +
                  0.0751 * tp3(ln_ggt_core, -1.379386) -
                  0.1235 * tp3(ln_ggt_core, -0.9416092) +
                  0.0540 * tp3(ln_ggt_core, 0.1655144) +
                  0.7062 * ln_ast_core -
                  0.2208 * tp3(ln_ast_core, -1.560648) +
                  0.7297 * tp3(ln_ast_core, -1.203973) -
                  0.9805 * tp3(ln_ast_core, -0.9942523) +
                  0.4716 * tp3(ln_ast_core, -0.4654975) -
                  0.6343 * ln_alt_core +
                  0.0257 * tp3(ln_alt_core, -1.89712) -
                  0.1750 * tp3(ln_alt_core, -1.272966) +
                  0.2565 * tp3(ln_alt_core, -0.836248) -
                  0.1072 * tp3(ln_alt_core, 0.008983292) +
                  male_core * (
                    -0.9480 +
                      0.4667 * ln_ggt_core +
                      0.0273 * tp3(ln_ggt_core, -1.966113) -
                      0.2593 * tp3(ln_ggt_core, -1.379386) +
                      0.3749 * tp3(ln_ggt_core, -0.9416092) -
                      0.1437 * tp3(ln_ggt_core, 0.1655144) -
                      0.1264 * ln_ast_core +
                      0.0890 * tp3(ln_ast_core, -1.560648) -
                      0.3117 * tp3(ln_ast_core, -1.203973) +
                      0.4582 * tp3(ln_ast_core, -0.9942523) -
                      0.2356 * tp3(ln_ast_core, -0.4654975) +
                      0.0964 * ln_alt_core -
                      0.0364 * tp3(ln_alt_core, -1.89712) +
                      0.1968 * tp3(ln_alt_core, -1.272966) -
                      0.2865 * tp3(ln_alt_core, -0.836248) +
                      0.1261 * tp3(ln_alt_core, 0.008983292)
                  ) +
                  age * (
                    0.0361 * ln_ggt_core +
                      0.0004 * tp3(ln_ggt_core, -1.966113) -
                      0.0028 * tp3(ln_ggt_core, -1.379386) +
                      0.0044 * tp3(ln_ggt_core, -0.9416092) -
                      0.0020 * tp3(ln_ggt_core, 0.1655144)
                  ) +
                  ln_ggt_core * (
                    -0.0320 * age +
                      0.0012 * tp3(age, 30.65) -
                      0.0039 * tp3(age, 39.35) +
                      0.0057 * tp3(age, 46.65) -
                      0.0030 * tp3(age, 56.18) +
                      0.8687 * ln_ast_core -
                      0.0083 * tp3(ln_ast_core, -1.560648) -
                      0.0658 * tp3(ln_ast_core, -1.203973) +
                      0.2139 * tp3(ln_ast_core, -0.9942523) -
                      0.1327 * tp3(ln_ast_core, -0.4654975) -
                      0.8478 * ln_alt_core -
                      0.2202 * tp3(ln_alt_core, -1.89712) +
                      1.1918 * tp3(ln_alt_core, -1.272966) -
                      1.7143 * tp3(ln_alt_core, -0.836248) +
                      0.7427 * tp3(ln_alt_core, 0.008983292)
                  ) +
                  ln_ast_core * (
                    -0.8745 * ln_ggt_core +
                      0.0089 * tp3(ln_ggt_core, -1.966113) +
                      0.0335 * tp3(ln_ggt_core, -1.379386) -
                      0.1793 * tp3(ln_ggt_core, -0.9416092) +
                      0.1275 * tp3(ln_ggt_core, 0.1655144) +
                      0.9496 * ln_alt_core +
                      0.1414 * tp3(ln_alt_core, -1.89712) -
                      0.7875 * tp3(ln_alt_core, -1.272966) +
                      1.1317 * tp3(ln_alt_core, -0.836248) -
                      0.4855 * tp3(ln_alt_core, 0.008983292)
                  ) +
                  ln_alt_core * (
                    0.8298 * ln_ggt_core +
                      0.1958 * tp3(ln_ggt_core, -1.966113) -
                      1.1645 * tp3(ln_ggt_core, -1.379386) +
                      1.6265 * tp3(ln_ggt_core, -0.9416092) -
                      0.6510 * tp3(ln_ggt_core, 0.1655144) -
                      0.9407 * ln_ast_core -
                      0.1428 * tp3(ln_ast_core, -1.560648) +
                      0.7842 * tp3(ln_ast_core, -1.203973) -
                      1.1137 * tp3(ln_ast_core, -0.9942523) +
                      0.4723 * tp3(ln_ast_core, -0.4654975)
                  )
)

master_nits <- master_nits %>%
  mutate(
    core_lp = lp_core,
    core_risk = 1 / (1 + exp(-core_lp)),
    core_percent = 100 * core_risk,
    core_z = as.numeric(scale(core_risk))
  )

n0 <- nrow(master_nits)
master_nits <- master_nits %>%
  filter(fib4 <= 10)
flow_counts <- bind_rows(flow_counts, add_flow("Final analytical cohort after exclusion of FIB-4 >10", n0, nrow(master_nits)))
