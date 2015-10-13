dir=rlx3

for mass in 0.50 0.60 0.70 0.80 0.90 1.00 1.10
do
    file="$dir"_s"$mass"/final.dat
    awk '{print $16}' $file \
        | sort -g | tail -n1 | awk '{print $1, mass}' mass=$mass
done

