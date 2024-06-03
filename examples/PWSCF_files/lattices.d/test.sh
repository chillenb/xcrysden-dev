PWI=../../../scripts/pwi2xsf.sh
PWO=../../../scripts/pwo2xsf.sh

pwis="
ibrav0_A.in
ibrav0_angs.in
ibrav0_bohr.in
ibrav0_celldm.in
ibrav1_celldm.in
ibrav2_A.in
ibrav2_A_crystal.in
ibrav2_celldm.in
ibrav2_celldm_crystal.in
ibrav3_A.in
ibrav3_A_angstrom.in
ibrav3_celldm.in
ibrav-3_celldm.in
ibrav4_celldm.in
ibrav5_celldm.in
ibrav-5_celldm.in
ibrav6_celldm.in
ibrav7_celldm.in
ibrav8_celldm.in
ibrav9_celldm.in
ibrav-9_celldm.in
ibrav91_celldm.in
ibrav10_celldm.in
ibrav11_celldm.in
ibrav12_celldm.in
ibrav-12_celldm.in
ibrav13_celldm.in
ibrav-13_celldm.in
ibrav14_abc.in
ibrav14_celldm.in
"

for pwi in $pwis
do

    pwo=${pwi%.in}.out
    
    echo "
------------------------------------------------------------------------
* pw.x input file = $pwi
------------------------------------------------------------------------
"
    pw.x < $pwi > $pwo

    $PWI $pwi > $pwi.xsf
    $PWO -ic $pwo > $pwo.xsf

    tkdiff $pwi.xsf $pwo.xsf
done

rm -f pwi2xsf.xsf_out *.old *.new
