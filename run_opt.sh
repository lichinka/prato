#
# Sync cubo
#
#rsync -avz ${HOME}/etc/dr/tun_par/prato/ luka@192.168.1.102:${HOME}/etc/dr/tun_par/prato/

#
# Start MPI jobs
#
mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage -pt ini_file=./parameters_gsm.ini tx_ini_sections=KPODVI1 : -np 1 --hostfile hostfile.local ./run_worker.sh
#mpirun --mca btl_tcp_if_include 192.168.1.0/24 -np 1 -host k100 r.coverage -p ini_file=./parameters_gsm.ini tx_ini_sections=KPODVI1 : -np 1 --hostfile hostfile.local ./run_worker.sh


############
## GPU-CPU comparison test-case follows

#
# Start MPI jobs for GPU
#
#mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage -g ini_file=./parameters.ini tx_ini_sections=CVELEA : -np 1 --hostfile hostfile.local ./run_worker.sh
#
# Spare the output of the GPU job
#
#cat /tmp/worker.log > /tmp/worker.gpu

#
# Start job for CPU
#
#r.ericsson inputDEM=dem100i_mo@PERMANENT clutter=clut100_mo@PERMANENT output=temp1 A0=38 A1=32 A2=-12 A3=0.1 coordinate=509535,135169 ant_height=30.65 radius=8 frequency=2040 --overwrite

#r.mapcalc diff=temp-temp1
#r.info diff
#r.colors map=diff color=elevation

#
# Spare the output of the CPU job
#
#cat /tmp/worker.log > /tmp/worker.cpu

#echo "***"
#echo "* $( wc -l /tmp/worker.cpu )"
#echo "* $( wc -l /tmp/worker.gpu )"
#echo "* GPU and CPU outputs differ in $( diff /tmp/worker.cpu /tmp/worker.gpu | wc -l )/$( wc -l /tmp/worker.cpu | cut -d' ' -f1 ) lines."
#echo "***"


