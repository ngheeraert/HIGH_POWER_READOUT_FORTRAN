nl=5
qb_ini=2
#cav_ini=1
cav_ini=0
ncs_ini=2
ncs_max=15
dvice=TRSM2
wc=7.415
g_qc=-0.10725
#g_qc=0.0000001
p0=0.000005
ad=0.006
#ad=0.00848
#ad=0
gama=0.0055585 #--  to get 1.2MHz decay
nbm=300
tmax=251 #-- rebound at t=400
offset_abs=1
offset_ang=0.5
p2_threshold=015
calc_error=1
merr=0.02
method=RK45
ar=0.0
lsfe=14
gcut=0.000001

n_dtadd_change=4
dtadd=2
dt_add_ini=$dtadd

dt=0.01
PD=0.03
ECD=0.03
median_pts=20

wd=7.419
#wd=7.415

bc_lf=2.4
bc_qb=4.461

bc_cav=7.419
bc_hf=12.6

bw=0.55
mr_lf=0.0
mr_qb=0.8
mr_hf=0.0


#merr_max=0.0001
#merr_int=0.0004

merr_array=(0.05 0.04 0.03 0.02 0.01)
#dtadd_array=(2)

errlim=0.0001
paramchar=mpol_nl$nl\_n$ncs_ini\_$ncs_max\_ad$Ad\_dtadd$dt_add\_mr$mr\_t$tmax\_errlim$errlim\_abs$offset_abs

#for dtadd in ${dtadd_array[@]}; 
for merr in ${merr_array[@]}; 
do

  	echo i=$i, ncs_max=$ncs_max, merr=$merr, dtadd=$dtadd, dt_add_ini=$dt_add_ini
  	newjobfile=$ncs_max\_$merr\_$dtadd
  	cp jobfile $newjobfile

	echo "./mpol -qb_ini $qb_ini -cav_ini $cav_ini -nl $nl -g_qc $g_qc -ncs_ini $ncs_ini -ncs_max $ncs_max -dvice $dvice -g_qc $g_qc -wc $wc -wd $wd -bc_lf $bc_lf -bc_qb $bc_qb -bc_cav $bc_cav -bc_hf $bc_hf -mr_lf $mr_lf -mr_qb $mr_qb -mr_hf $mr_hf -p0 $p0 -ad $ad -bw $bw -dt $dt -dt_add $dtadd -gama $gama -nbm $nbm -tmax $tmax -errlim $errlim -error_thr $merr -offset_abs $offset_abs -offset_ang $offset_ang -p2_threshold $p2_threshold -logfile 2 -print_delay $PD -error_calc_delay $ECD -calc_error $calc_error -method $method -adding_ratio $ar -job_nb -1 -median_pts $median_pts -n_dtadd_change $n_dtadd_change -dt_add_ini $dt_add_ini -lsfe $lsfe -gcut $gcut" >> $newjobfile

	echo "mv ../MPOL_\$tpdir/data/* \$PBS_O_WORKDIR/data/" >> $newjobfile
	echo "rm -r ../MPOL_\$tpdir" >> $newjobfile

	qsub $newjobfile
	rm $newjobfile

done
