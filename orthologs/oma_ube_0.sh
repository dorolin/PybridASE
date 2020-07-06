sbatch --nodes=1 --mail-user=me@work --mail-type=end,fail --job-name="Oma0" --output=%x-%j.out --time=01:30:00 --mem-per-cpu=5000M <<EOF
#!/bin/sh
echo "Starting at `date`"
OMA -c
EOF
