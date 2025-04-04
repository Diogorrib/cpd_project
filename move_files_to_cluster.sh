IST_ID="ist1102484"
CLUSTER="cluster.rnl.tecnico.ulisboa.pt"
CLUSTER_DIR="/mnt/cirrus/users/8/4/$IST_ID/CPD"

# Copy the files to the cluster
scp -r test_mpi.sh collect_results_and_time.sh mpi samples $IST_ID@$CLUSTER:$CLUSTER_DIR
