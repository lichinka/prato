#
# Sync cubo
#
rsync -avz /home/luka/etc/dr/tun_par/prato/ luka@192.168.1.101:/home/luka/etc/dr/tun_par/prato/

#
# Start MPI jobs
#
mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage -t ini_file=./parameters.ini tx_ini_sections=CVELEA : -np 1 --hostfile hostfile.local ./run_worker.sh
