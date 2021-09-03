# sleep_hormone_MR

Two-sample Mendelian randomization to investigate causal relationships between sleep traits, hormone traits and risk of breast/prostate cancer. 

# R scripts 
> uvMR_breast_and_prostate_cancer.R

This script can be edited and used to run:
Chronotype > Breast Cancer (overall)
Total Testosterone > Breast Cancer (overall)
Bioavailable Testosterone > Breast Cancer (overall)
SHBG > Breast Cancer (overall)
Oestradiol > Breast Cancer (overall)

Chronotype > Prostate Cancer
Total Testosterone > Prostate Cancer
Bioavailable Testosterone > Prostate Cancer
SHBG > Prostate Cancer
Oestradiol > Prostate Cancer

> bdMR_sleep_and_hormone_traits.R 

This script can be edited and used to run:
Chronotype <> Total Testosterone (female)
Chronotype <> Bioavailable Testosterone (female)

Chronotype <> Bioavailable Testosterone (male)

Both exposure and outcome data in bdMR are from UK Biobank, as such, Steiger filtering is applied as standard.

> mvMR_breast_and_prostate_cancer.R

This script can be edited and used to run:
Chronotype + Total Testosterone > Breast Cancer
Chronotype + Bioavailable Testosterone > Breast Cancer

Chronotype + Bioavailable Testosterone > Prostate Cancer

