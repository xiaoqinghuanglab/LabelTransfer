
library(tidyverse)
library(ggpubr)
library(rstatix)

rosmap_blood_clinical = read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/ANMerge/AD_LabelTransfer_GitHub/Data/ROSMAP_blood_clinical_monocyte_ACTL_04172023.csv")
rosmap_blood_ge = read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/ANMerge/AD_LabelTransfer_GitHub/Data/ROSMAP_blood_gene_expression_monocyte_ACTL_04172023.csv")
rosmap_blood_clinical$AD_Subtype = factor(rosmap_blood_clinical$Group1, levels = c("Asym AD", "Control", "Typical AD", "Low-NFT AD"))

# Statistical test
stat.test <- rosmap_blood_clinical %>% t_test(cts_mmse30_lv ~ AD_Subtype)
stat.test
bxp <- ggboxplot(rosmap_blood_clinical, x = "AD_Subtype", y = "cts_mmse30_lv", color = "AD_Subtype", ylim=c(0, 52)) 
# Box plot
stat.test <- stat.test %>% add_xy_position(x = "AD_Subtype")
bxp +
  ggtitle("ROSMAP blood MMSE by subtypes") +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label = p.adj.signif), # label = signif(p, 2)),
    data = stat.test, y.position = 35, step.increase = 0.1
  )

## ADNI
adni_blood_clinical = read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/Github/Data/ADNI_NewModel_NewLabel_Clinical_Info.csv")
adni_blood_clinical = adni_blood_clinical[!adni_blood_clinical$NewLabel == "Others", ]
adni_blood_clinical$NewLabel = factor(adni_blood_clinical$NewLabel, levels = c("Asym AD", "Control", "Typical AD", "Low-NFT AD"))

# Statistical test
stat.test <- adni_blood_clinical %>% t_test(MMSE_Close ~ NewLabel)
stat.test
bxp <- ggboxplot(adni_blood_clinical, x = "NewLabel", y = "MMSE_Close", color = "NewLabel", ylim=c(0, 52)) 
# Box plot
stat.test <- stat.test %>% add_xy_position(x = "NewLabel")
bxp +
  ggtitle("ADNI blood MMSE by subtypes") +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label = p.adj.signif), # label = signif(p, 2)),
    data = stat.test, y.position = 35, step.increase = 0.1
  )

## ANMerge
anmerge_blood_clinical = read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/Github/Data/anmerge_NewModel_NewLabel_Clinical_Info.csv")
anmerge_blood_clinical = anmerge_blood_clinical[!anmerge_blood_clinical$NewLabel == "Others", ]
anmerge_blood_clinical$NewLabel = factor(anmerge_blood_clinical$NewLabel, levels = c("Asym AD", "Control", "Typical AD", "Low-NFT AD"))
# Statistical test
# anmerge_blood_clinical = anmerge_blood_clinical[anmerge_blood_clinical$MMSE != Inf, ]
stat.test <- anmerge_blood_clinical[anmerge_blood_clinical$MMSE != Inf, ] %>% t_test(MMSE ~ NewLabel)
stat.test
bxp <- ggboxplot(anmerge_blood_clinical, x = "NewLabel", y = "MMSE", color = "NewLabel", ylim=c(0, 52)) 
# Box plot
stat.test <- stat.test %>% add_xy_position(x = "NewLabel")
bxp +
  ggtitle("ANMerge blood MMSE by subtypes") +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label = p.adj.signif), # label = signif(p, 2)),
    data = stat.test, y.position = 35, step.increase = 0.1
  )


######## old label
## ADNI
adni_blood_clinical = read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/Github/Data/ADNI_NewModel_NewLabel_Clinical_Info.csv")
#adni_blood_clinical = adni_blood_clinical[!adni_blood_clinical$NewLabel == "Others", ]
adni_blood_clinical$Label = factor(adni_blood_clinical$Label, levels = c("Asym AD", "Control", "Typical AD", "Low-NFT AD"))

# Statistical test
stat.test <- adni_blood_clinical %>% t_test(MMSE_Close ~ Label)
stat.test
bxp <- ggboxplot(adni_blood_clinical, x = "Label", y = "MMSE_Close", color = "Label", ylim=c(0, 52)) 
# Box plot
stat.test <- stat.test %>% add_xy_position(x = "Label")
bxp +
  ggtitle("ADNI blood MMSE by subtypes: OLD Label") +
  geom_bracket(
    aes(xmin = group1, xmax = group2, label = p.adj.signif), # label = signif(p, 2)),
    data = stat.test, y.position = 35, step.increase = 0.1
  )

