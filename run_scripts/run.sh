
snakemake --configfile config/config_t2t_v2.0.yaml --cores $(nproc) --rerun-triggers mtime \
  $@

exit

snakemake --configfile config/config.yaml --cores $(nproc) --rerun-triggers mtime \
  $@
