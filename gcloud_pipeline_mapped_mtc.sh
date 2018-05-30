#!/bin/bash

my_bucket="gbsc-gcp-lab-snyder-user-ahorning"
my_project="gbsc-gcp-lab-snyder" ###This is your project name that contains the BigQuery dataset.
# my-bigquery-dataset: FAP_HaplotypeCaller_gVCFs ###This is the dataset ID you created in the prior step.
# my-bigquery-table: gVCFs_table This can be any ID you like such as “platinum_genomes_variants”.

gcloud config set project gbsc-gcp-lab-snyder

gcloud alpha genomics pipelines run \
    --project "${my_project}" \
    --pipeline-file vcf_to_bigquery_mapped_mtc.yaml \
    --logging gs://"${my_bucket}"/temp/runner_logs \
    --zones us-west1-b \
    --service-account-scopes https://www.googleapis.com/auth/bigquery


# Updated property [core/project].
# Running [operations/EL-L_q-iLBjD_vCsitfd1yUgh73Nq78ZKg9wcm9kdWN0aW9uUXVldWU].
#gcloud alpha genomics operations describe EL-L_q-iLBjD_vCsitfd1yUgh73Nq78ZKg9wcm9kdWN0aW9uUXVldWU