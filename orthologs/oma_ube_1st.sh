sbatch --array=1-100 --nodes=1 --mail-user=me@work --mail-type=end,fail --job-name="Oma1st" --output=%x-%j.out --time=3-23:30:00 --mem-per-cpu=2150M <<EOF
#!/bin/sh
export NR_PROCESSES=100
echo "Starting at `date`"
OMA -s -W 342000
EOF

## grep "success" in log files to check for completion
