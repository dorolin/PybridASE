sbatch --nodes=1 --mail-user=me@work --mail-type=end,fail --job-name="Oma2nd" --output=%x-%j.out --time=23:55:00 --mem-per-cpu=15G <<EOF
#!/bin/sh
echo "Starting at `date`"
OMA
EOF

## look for "syntax error" and "unexpected end of file" in log files: 
##   in that case, delete all files in Cache/AllAll/ for species pair 
##   that comes after the error and re-run 'oma_ube_1st.ch' and
##   'oma_ube_2nd.sh' and check that re-run removes error
