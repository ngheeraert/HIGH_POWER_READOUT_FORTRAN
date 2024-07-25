nl=2
qb_ini=2
cav_ini=0
ncs_ini=2
ncs_max=20
dvice=TRSM3
wc=7.5
#g_qc=0
g_qc=0.25
p0=0.000005
ad=0.0
#ad=0.05
#ad=0.2
#ad=0.4
#alpha=0.00
alpha=0.001071759819112394
#alpha_LF=0.0
alpha_LF=0.000025
cutoff_LF=6.5
nbm=500
#nbm=1000
tmax=25 #-- rebound at t=120
offset_abs=1.5
offset_ang=1
p2_threshold=015	
calc_error=1
#merr=0.0005
method=RK45
ar=0.0
lsfe=14
ncut=0.000001

n_dtadd_change=2
#dtadd=0.2
#dt_add_ini=0.3
dtadd=0.2
dt_add_ini=0.02

dt=0.005
#PD=0.025
PD=0.005
ECD=0.005
median_pts=20

wd=7.485

bc_lf=2.3
bc_qb=5.2
bc_cav=6
bc_hf=12.6

#bw=4
#bw=8
#bw=32
bath_type='O'
bw=12
mr_lf=0.0
mr_qb=0.0
mr_hf=0.0

#merr_max=0.0001
#merr_int=0.0004

merr_array=(0.01)
#dtadd_array=(2)

errlim=0.00001
paramchar=mpol_nl$nl\_n$ncs_ini\_$ncs_max\_ad$Ad\_dtadd$dt_add\_mr$mr\_t$tmax\_errlim$errlim\_abs$offset_abs

#for dtadd in ${dtadd_array[@]}; 
for merr in ${merr_array[@]}; 
do

  	echo i=$i, ncs_max=$ncs_max, merr=$merr, dtadd=$dtadd, dt_add_ini=$dt_add_ini

  	newjobfile=$ncs_max\_$merr\_$dtadd
  	cp jobfile $newjobfile


	echo "./mpol -qb_ini $qb_ini -cav_ini $cav_ini -nl $nl -ncs_ini $ncs_ini -ncs_max $ncs_max -dvice $dvice -g_qc $g_qc -wc $wc -wd $wd -bc_lf $bc_lf -bc_qb $bc_qb -bc_cav $bc_cav -bc_hf $bc_hf -mr_lf $mr_lf -mr_qb $mr_qb -mr_hf $mr_hf -p0 $p0 -ad $ad -bw $bw -dt $dt -dt_add $dtadd -al $alpha -al_LF $alpha_LF -cutoff_LF $cutoff_LF -nbm $nbm -tmax $tmax -errlim $errlim -error_thr $merr -offset_abs $offset_abs -offset_ang $offset_ang -p2_threshold $p2_threshold -logfile 2 -print_delay $PD -error_calc_delay $ECD -calc_error $calc_error -adding_ratio $ar -median_pts $median_pts -n_dtadd_change $n_dtadd_change -dt_add_ini $dt_add_ini -lsfe $lsfe -ncut $ncut -bath_type $bath_type"  >> $newjobfile


	echo "mv ../MPOL_\$tpdir/data/* \$PBS_O_WORKDIR/data/" >> $newjobfile
	echo "rm -r ../MPOL_\$tpdir" >> $newjobfile

	qsub $newjobfile
	rm $newjobfile

done

