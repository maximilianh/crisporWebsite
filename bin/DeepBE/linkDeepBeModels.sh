#!/bin/bash

# This script removes the empty model directories and files in DeepBE,
# and replace them with symlinks from DeepBE/models.
# Run this script after downloadDeepBeModels.sh

baseDir="$(dirname "$(realpath "$0")")"

# array mapping models to their destination (there should be a cleaner way to handle this)
declare -A mapDirs=(
[DeepBE_17m_model]="DeepBE/ABE8.17m+V106W"
[DeepBE_8e_model]="DeepBE/ABE8e(V106W)"
[DeepBE_Bi_model]="DeepBE/APOBEC-nCas9-Ung"
[DeepBE_CGBE1_model]="DeepBE/CGBE1"
[DeepBE_mini_model]="DeepBE/miniCGBE1"
[DeepBE_Ss_model]="DeepBE/SsAPOBEC3B"
[DeepBE_YE1_model]="DeepBE/YE1-BE4max"
[DeepNG-BE_17m_model]="DeepNG-BE/SpCas9-NG-ABE8.17m+V106W"
[DeepNG-BE_8e_model]="DeepNG-BE/SpCas9-NG-ABE8e(V106W)"
[DeepNG-BE_Bi_model]="DeepNG-BE/SpCas9-NG-APOBEC-nCas9-Ung"
[DeepNG-BE_CGBE1_model]="DeepNG-BE/SpCas9-NG-CGBE1"
[DeepNG-BE_mini_model]="DeepNG-BE/SpCas9-NG-miniCGBE1"
[DeepNG-BE_Ss_model]="DeepNG-BE/SpCas9-NG-SsAPOBEC3B"
[DeepNG-BE_YE1_model]="DeepNG-BE/SpCas9-NG-YE1-BE4max"
[PAM_variant_NG_model.h5]="PAM/PAM_variant_NG"
[PAM_variant_NRCH_model.h5]="PAM/PAM_variant_NRCH"
[PAM_variant_NRRH_model.h5]="PAM/PAM_variant_NRRH"
[PAM_variant_NRTH_model.h5]="PAM/PAM_variant_NRTH"
[PAM_variant_Sc++_model.h5]="PAM/PAM_variant_Sc++"
[PAM_variant_SpCas9_model.h5]="PAM/PAM_variant_SpCas9"
[PAM_variant_SpG_model.h5]="PAM/PAM_variant_SpG"
[PAM_variant_SpRY_model.h5]="PAM/PAM_variant_SpRY"
[PAM_variant_VRQR_model.h5]="PAM/PAM_variant_VRQR"
)

# remove placehodler model files and dirs from the repo
find "$baseDir"/PAM -type f -name "*.h5" -exec rm {} ";"
find "$baseDir"/DeepBE -type d -name "*model" -exec rm -r {} ";"
find "$baseDir"/DeepNG-BE -type d -name "*model" -exec rm -r {} ";"

# replace them with symlinks
for dir in "$baseDir"/models/*; do
    [ -d "$dir" ] || continue
    for model in $dir/*; do
        modelName="$(basename $model)"
        destDir=$baseDir/"${mapDirs[$modelName]}"
        ln -sf "$model" "$destDir"/"$modelName" 
    done
done

# symlink the models in models/DeepBE/PAM to DeepBE/DeepBE/
for dir in "$baseDir"/DeepBE/*; do
    [ -d "$dir" ] || continue
    for modelDir in "$baseDir"/models/PAM/*; do
        PamModel=$(find $modelDir -type f -name "*.h5")
        PamModelName="$(basename "$PamModel")"
        ln -sf "$PamModel" "$dir/$PamModelName"
    done
done
