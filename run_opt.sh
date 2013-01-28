#
# Sync cubo
#
#rsync -avz ${HOME}/etc/dr/tun_par/prato/ luka@192.168.1.100:${HOME}/etc/dr/tun_par/prato/

#
# Start MPI jobs
#
mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage -gt ini_file=./parameters.ini tx_ini_sections=CVELEA : -np 1 --hostfile hostfile.local ./run_worker.sh


############
## GPU-CPU comparison test-case follows

#
# Start MPI jobs
#
#mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage -g ini_file=./parameters.ini tx_ini_sections=CTOPOLA : -np 1 --hostfile hostfile.local ./run_worker.sh
#
# Spare the output of the GPU job
#
#cat /tmp/worker.log > /tmp/worker.gpu

#
# Start MPI jobs
#
#mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage -p ini_file=./parameters.ini tx_ini_sections=CTOPOLA : -np 1 --hostfile hostfile.local ./run_worker.sh
#
# Spare the output of the CPU job
#
#cat /tmp/worker.log > /tmp/worker.cpu

#echo "***"
#echo "* GPU and CPU outputs differ in $( diff /tmp/worker.cpu /tmp/worker.gpu | wc -l )/$( wc -l /tmp/worker.cpu | cut -d' ' -f1 ) lines."
#echo "***"


