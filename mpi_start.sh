#
# Sync cubo
#
rsync -avz /home/luka/etc/dr/tun_par/prato/ luka@192.168.1.102:/home/luka/etc/dr/tun_par/prato/

#
# Start MPI jobs
#
mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage -p ini_file=/home/luka/etc/dr/tun_par/prato/src/parameters.ini tx_ini_sections=CTOPOLA,CVELCEA,SBANOVA,CVELCEB,CTOPOLA,SBANOVA : -np 2 --hostfile hostfile.local ./run_worker.sh
