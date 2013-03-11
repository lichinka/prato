#
# Sync cubo
#
#rsync -avz ${HOME}/etc/dr/tun_par/prato/ luka@192.168.1.102:${HOME}/etc/dr/tun_par/prato/

#
# Start MPI jobs
#
#mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage $1 ini_file=./parameters_gsm.ini tx_ini_sections=KPODVI1,KPODVI1,KPODVI1 : -np 1 --hostfile hostfile.local ./run_worker.sh
mpirun --mca btl_tcp_if_include 192.168.1.0/24 --mca btl tcp,sm,self -np 1 -host localhost r.coverage $1 ini_file=./parameters_umts.ini tx_ini_sections=CSOSTA,CSOSTB,CSOSTC : -np 1 --hostfile hostfile.local ./run_worker.sh


############
## GPU-CPU comparison test-case follows

#
# Start MPI jobs for GPU
#
#mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage -g ini_file=./parameters_gsm.ini tx_ini_sections=KPODVI1 : -np 1 --hostfile hostfile.local ./run_worker.sh
#
# Spare the output of the GPU job
#
#mv /tmp/worker.log /tmp/worker.gpu

#
# Start job for CPU
#
#mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage -p ini_file=./parameters_gsm.ini tx_ini_sections=KPODVI1 : -np 1 --hostfile hostfile.local ./run_worker.sh

#
# Spare the output of the CPU job
#
#mv /tmp/worker.log /tmp/worker.cpu

#echo "***"
#echo "* $( wc -l /tmp/worker.cpu )"
#echo "* $( wc -l /tmp/worker.gpu )"
#echo "* GPU and CPU outputs differ in $( diff /tmp/worker.cpu /tmp/worker.gpu | wc -l )/$( wc -l /tmp/worker.cpu | cut -d' ' -f1 ) lines."
#echo "***"

#r.mapcalc diff=temp-temp1
#r.info diff
#r.colors map=diff color=elevation

