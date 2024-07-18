nl=5
qb_ini=1
#cav_ini=1
cav_ini=0
ncs_ini=2
ncs_max=2
dvice=QUTR2
wc=7.415
#g_qc_TRANS=-0.10725
#g_qc=0.000001
g_qc=0.00369
p0=0.000005
#ad=0.003
ad=0.009
#ad=0.0
gama=0.0055585 #--  to get 1.2MHz decay
nbm=300
tmax=200 #-- rebound at t=400
offset_abs=1
offset_ang=0.3
p2_threshold=015
calc_error=0
merr=-0.001
method=RK45
ar=0.0
lsfe=16
gcut=0.000001

n_dtadd_change=4
dtadd=50
dt_add_ini=$dtadd

dt=0.01
PD=0.03
ECD=0.03
median_pts=20

wd=7.4211
#wd=7.415

bc_lf=2.4
bc_qb=4.461

bc_cav=7.4211
bc_hf=12.6

#bw=0.55
bw=2.2
mr_lf=0.0
mr_qb=0.8
mr_hf=0.0

errlim=0.000000001
paramchar=mpol_nl$nl\_n$ncs_ini\_$ncs_max\_ad$Ad\_dtadd$dt_add\_mr$mr\_t$tmax\_errlim$errlim\_abs$offset_abs

echo i=$i, ncs_max=$ncs_max, merr=$merr, dtadd=$dtadd, dt_add_ini=$dt_add_ini
#newjobfile=$ncs_max\_$merr\_$dtadd
#cp jobfile $newjobfile

echo "./mpol -qb_ini $qb_ini -cav_ini $cav_ini -nl $nl -g_qc $g_qc -ncs_ini $ncs_ini -ncs_max $ncs_max -dvice $dvice -g_qc $g_qc -wc $wc -wd $wd -bc_lf $bc_lf -bc_qb $bc_qb -bc_cav $bc_cav -bc_hf $bc_hf -mr_lf $mr_lf -mr_qb $mr_qb -mr_hf $mr_hf -p0 $p0 -ad $ad -bw $bw -dt $dt -dt_add $dtadd -gama $gama -nbm $nbm -tmax $tmax -errlim $errlim -error_thr $merr -offset_abs $offset_abs -offset_ang $offset_ang -p2_threshold $p2_threshold -logfile 2 -print_delay $PD -error_calc_delay $ECD -calc_error $calc_error -method $method -adding_ratio $ar -job_nb -1 -median_pts $median_pts -n_dtadd_change $n_dtadd_change -dt_add_ini $dt_add_ini -lsfe $lsfe -gcut $gcut" #>> $newjobfile

#echo "mv ../MPOL_\$tpdir/data/* \$PBS_O_WORKDIR/data/" >> $newjobfile
#echo "rm -r ../MPOL_\$tpdir" >> $newjobfile

#qsub $newjobfile
#rm $newjobfile

