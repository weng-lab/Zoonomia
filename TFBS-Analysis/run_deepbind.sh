#!/bin/bash

acc=$1


module load gcc/4.8.1
module load cuda/7.0.28
source activate deepbind 

export LD_LIBRARY_PATH=/project/umw_zhiping_weng/andrewsg/gifford/deepbind-docker/code/libs/deepity/build/release/bin/:/project/umw_zhiping_weng/andrewsg/gifford/deepbind-docker/code/libs/smat/build/release/bin/:/project/umw_zhiping_weng/andrewsg/gifford/deepbind-docker/code/libs/kangaroo/build/release/bin:$LD_LIBRARY_PATH

cd /project/umw_zhiping_weng/andrewsg/gifford/deepbind-docker/code
## Run deepbind on ENCDOE ChIP-seq data
#python deepbind_train_encode.py all calib,train,test,report -i ../data/motif_discovery/$acc/ -o ../data/motif_discovery/$acc/deepbind
#python deepbind_train_encode.py all calib -i ../data/motif_discovery/$acc/ -o ../data/motif_discovery/$acc/deepbind
#python deepbind_train_encode.py all train -i ../data/motif_discovery/$acc/ -o ../data/motif_discovery/$acc/deepbind
#python deepbind_train_encode.py all test -i ../data/motif_discovery/$acc/ -o ../data/motif_discovery/$acc/deepbind
#python deepbind_train_encode.py all report -i ../data/motif_discovery/$acc/ -o ../data/motif_discovery/$acc/deepbind

#python deepbind_train_encode.py all calib -i ../data/encode/$acc/ -o ../data/encode/$acc/deepbind
#python deepbind_train_encode.py all train -i ../data/encode/$acc/ -o ../data/encode/$acc/deepbind
python deepbind_train_encode.py all test -i ../data/encode/$acc/ -o ../data/encode/$acc/deepbind
#python deepbind_train_encode.py all report -i ../data/encode/$acc/ -o ../data/encode/$acc/deepbind


## Run DeepBind on HT-SELEX
#python deepbind_train_selex.py --jolma calib,train,test,report -i ../data/selex/jolma/$acc -o ../data/selex/jolma/$acc/deepbind_out
#python deepbind_train_selex.py --jolma test,report -i ../data/selex/jolma/$acc -o ../data/selex/jolma/$acc/deepbind_out
source deactivate 