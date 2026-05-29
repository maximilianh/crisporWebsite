#!/bin/bash

# this script downlads the models for DeepBE
# see https://github.com/NahyeKim/DeepBE/releases/tag/version1
# run it once on new installations

baseDir="$(dirname "$(realpath "$0")")"
models=(a b c d e f g h i j)

mkdir -p "$baseDir"/models

for model in "${models[@]}"; do
echo "-------------- downloading part "$model" ---------------"
    curl --progress-bar --output-dir "$baseDir"/models -O -L https://github.com/NahyeKim/DeepBE/releases/download/version1/model_v1"$model"
done

echo "---------- merging and extracting data ----------"
modelArchive="$baseDir"/models/model_v1.tar.gz
cat "$baseDir"/models/model_v1* > "$modelArchive"
tar xvzf "$modelArchive" --no-same-owner -C "$baseDir"/models

echo "-------------- cleaning directory ---------------"
rm "$modelArchive"
for model in "${models[@]}"; do
    rm "$baseDir"/models/model_v1"$model"
done
