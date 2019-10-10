import dalmatian as dm
import pandas as pd
import Helper
import os
import numpy as np
import ipdb
import TerraFunction as terra


def getCNHeatMap(workspace1="CCLF_TSCA_2_0_2", namespace1="nci-mimoun-bi-org",
                 list_sample_ids = None,
                 samplesetname = None,
                 outputloc= None
                 ):
  #####
  ### NEED TO COMPLETE / FIX ###
  #####
  # produce a copy number heatmap for all the samples in list_sample_ids
  wm1 = dm.WorkspaceManager(namespace1, workspace1)
  wm1.update_sample_set(samplesetname, list_sample_ids) # create new sample set containing the samples in the list
  pathto_cnvpng = "cnv_calls_tn_img"
  # run the __ workflow to produce a merged, labeled copy number heat map
  # need to get access permissions to run anything in TSCA workspace...
  # submission_id = wm1.create_submission(CGAmutations['name'], samplesetname,'sample_set',expression='this.samples')
  # terra.waitForSubmission(wm1, submission_id)

  # copy the CN plot to the output directory
  sample_set = wml.get_sample_sets()
  cnfile = sample_set["sample_set_id"==samplesetname][pathto_cnvpng]
  os.system('gsutil cp ' + cnfile + ' ' + outputloc + 'copy_number_heatmap.png')



def getReport(workspace1="CCLF_Targeted", namespace1="nci-mimoun-bi-org", # CCLF_Targeted instead of CCLF_TSCA_2_0_2
              pathto_cnvpng='segmented_copy_ratio_img',
              pathto_stats='sample_statistics',
              pathto_snv='filtered_variants',
              pathto_seg='cnv_calls',
              workspacewes='CCLF_WES', namespacewes='nci-mimoun-bi-org',
              pathto_cnvpng_wes='segmented_copy_ratio_img',
              pathto_stats_wes='sample_statistics',
              pathto_snv_wes='mafLite',
              pathto_seg_wes='tumor_seg',
              is_from_pairs=True,
              datadir='gs://cclf_results/targeted/kim_sept/',
              tempdir='temp/cclfmerge/',
              specificlist=None,
              specificlist_disease=None
              ):
  """
  """
  # specificlist = filename/filepath of file containing a list of participant IDs
  # if we get many condition: merge and display the one that has been selected
  ## Neekesh wants all of them to be included, but to add additional mark to the one that has been selected

  # if we get WES data: add on the side as additional data and merge CNplot
  # if we get primary: add on the side as additional data and merge CNplot
  # if have paired data, include normal in CNplot
  print('you need to be on macOS for now')
  wm1 = dm.WorkspaceManager(namespace1, workspace1)
  wm_wes = dm.WorkspaceManager(namespacewes, workspacewes)
  if type(specificlist) is str:
    # we consider it is a filename
    specificlist = pd.read_csv(specificlist).tolist()
  elif specificlist is None:
    print("taking all samples from the workspace")
    specificlist = wm1.get_sample_sets().index

  sample = wm1.get_samples()
  sample_wes = wm_wes.get_samples()
  if is_from_pairs:
    pair = wm1.get_pairs()
    pair_wes = wm_wes.get_pairs()
  sample_part = sample.participant.tolist()
  samplewes_part = sample_wes.participant.tolist()

  # get data for each participant ID (e.g. PEDS182)
  for val in specificlist:
    found = False
    found_TSCA = False
    mutfile = pd.DataFrame()
    cnfile = pd.DataFrame()
    images = [] # store all the CN images for a participant
    ext_ids = [] # will collect all the external_ids we have for a participant

    # get data from targeted probes (TSCA, TWIST)
    if val in sample_part:
      sample_subset = sample[sample.participant == val]
      samples_for_CN_heat_TSCA = [] # sample set to make CN heat map for; want list of sample_id
      for k, condition in sample_subset.iterrows():
        if condition['sample_type'] == "Tumor":
          found = True
          found_TSCA = True
          samples_for_CN_heat_TSCA += [k] # k is the sample_id
          external_id = condition['external_id_validation']
          ext_ids += [external_id]
          cond_name = condition['media'] # problem w/ cond_name: can have multiple samples with the same participant_id and media.
          outputloc = datadir + condition['primary_disease'].replace(' ', '_') + '/' + val + '/'
          os.system('gsutil cp ' + condition[pathto_seg] + ' ' + tempdir + 'copy_number.tsv')
          os.system('gsutil cp ' + condition[pathto_cnvpng] + ' ' + tempdir + external_id + '_copy_number_map.png')
          # annotate image with sample information
          imagedir = tempdir + external_id + '_copy_number_map.png'
          text = condition["external_id_validation"].replace("_", " ") + " (TSCA)"
          Helper.addTextToImage(imagedir, text, outputpath=imagedir)
          images.append(tempdir + external_id + '_copy_number_map.png')
          os.system('gsutil cp ' + imagedir + ' ' + outputloc + external_id + '_copy_number_map.png')

          os.system('gsutil cp ' + condition[pathto_stats] + ' ' + outputloc + external_id + '_sample_statistics.txt')
          if is_from_pairs:
            ipdb.set_trace()
            pair_subset = pair[pair["case_sample"] == k]
            snvs = pair_subset[pathto_snv]
            for snv in snvs:
              if snv is not np.nan:
                os.system('gsutil cp ' + snv + ' ' + tempdir + 'mutations.tsv')
                # add first matched normal to list if exists (aka add the normal's sample_id)
                try:
                    samples_for_CN_heat_TSCA += [next(id for id in pair_subset["control_sample"] if id not in ["NA",np.nan])]
                except: # only NA or nans; no real matched normal
                    continue
                break
          else:
            os.system('gsutil cp ' + condition[pathto_snv_unpaired] + ' ' + tempdir + 'mutations.tsv')

          # add culture condition and external_id columns to the tables
          mut = pd.read_csv(tempdir + 'mutations.tsv', sep='\t')
          mut['condition'] = cond_name
          mut['external_id'] = external_id
          mutfile = mutfile.append(mut)
          cn = pd.read_csv(tempdir + 'copy_number.tsv', sep='\t')

          cn['condition'] = cond_name
          cn['external_id'] = external_id
          # cn['condition'] = cond_JKBname
          cnfile = cnfile.append(cn)
      if found:
        cnfile.to_csv(tempdir + 'cn.tsv', sep='\t')
        mutfile.to_csv(tempdir + 'mut.tsv', sep='\t')
        os.system('gsutil cp ' + tempdir + 'cn.tsv ' + outputloc + 'copy_number.tsv')
        os.system('gsutil cp ' + tempdir + 'mut.tsv ' + outputloc + 'mutation.tsv')
        # ext_ids.to_csv(tempdir + 'all_external_ids.tsv', sep='\t')
        # os.system('gsutil cp ' + tempdir + 'all_external_ids.tsv ' + outputloc + 'all_external_ids.tsv')

      # get CN plot for matched normal if a matched normal exists for a given participant (assuming only 1 distinct matched normal exists)
      matched_normal_ids = pair[pair["case_sample"] == k]["control_sample"]
      for normal_id in matched_normal_ids:
        if normal_id is not np.nan:
          # get the image
          normal_sample = sample[sample.index == normal_id] # indexed by sample_id
          if normal_sample[pathto_cnvpng][normal_id] not in ["NA",np.nan]:
            os.system('gsutil cp ' + normal_sample[pathto_cnvpng][normal_id] + ' ' + tempdir + normal_id + '_copy_number_map.png')

            # annotate matched normal CN image
            imagedir = tempdir + normal_id + '_copy_number_map.png'
            if os.path.exists(imagedir):
              text = normal_sample["external_id_validation"][normal_id].replace("_", " ") + " (TSCA, matched normal)"
              Helper.addTextToImage(imagedir, text, outputpath=imagedir)
              images.append(tempdir + normal_id + '_copy_number_map.png')
              os.system('gsutil cp ' + imagedir + ' ' + outputloc + normal_id + '_copy_number_map.png')
            break
      # # annotate matched normal CN image
      # imagedir = tempdir + normal_id + '_copy_number_map.png'
      # if os.path.exists(imagedir):
      #   text = normal_sample["external_id_validation"][normal_id].replace("_", " ") + " (TSCA, matched normal)"
      #   Helper.addTextToImage(imagedir, text, outputpath=imagedir)
      #   images.append(tempdir + normal_id + '_copy_number_map.png')
      #   os.system('gsutil cp ' + imagedir + ' ' + outputloc + normal_id + '_copy_number_map.png')

    # get data from CCLF_WES
    if val in samplewes_part:
      sample_subset_wes = sample_wes[sample_wes.participant == val]
      samples_for_CN_heat_WES = [] # sample set to make CN heat map for; want list of sample_id
      for k, wes in sample_subset_wes.iterrows():
        if wes['sample_type'] == "Tumor":
          found = True
          samples_for_CN_heat_WES += [k] # k is the sample_id
          # get primary disease information
          if found_TSCA: # pull disease info from TSCA
            primary = condition['primary_disease'].replace(' ', '_') if 'primary_disease' in condition else 'unknown'
          elif type(specificlist_disease) is str:
            # we consider it is a filename; file with info including primary_disease and participant_id
            specificlist_disease_df = pd.read_csv(specificlist_disease, index_col="participant_id")
            # grab all the disease information for the participant
            pri_disease_info = specificlist_disease_df.loc[val, "primary_disease"]
            if type(pri_disease_info) is not str:
                primary = pri_disease_info.iloc[0].replace(' ', '_') if 'primary_disease' in specificlist_disease_df else 'unknown'
            else:
                primary = pri_disease_info.replace(' ', '_') if 'primary_disease' in specificlist_disease_df else 'unknown'
          # if primary_disease info is added to CCLF_WES sometime, then just use the following:
          # primary = wes['primary_disease'].replace(' ', '_') if 'primary_disease' in wes else 'unknown'

          outputloc = datadir + primary + '/' + val + '/'
          external_id = wes['external_id_capture']
          ext_ids += [external_id]

          # get copy number information and plot
          ## have to check for np.nan; this means that we have data but haven't run the pipeline yet!!
          # os.system('gsutil cp ' + wes[pathto_seg_wes] + ' ' + outputloc + 'wes_copy_number.tsv')
          # os.system('gsutil cp ' + wes[pathto_cnvpng_wes] + ' ' + tempdir + 'wes_copy_number_map.pdf')

          if wes[pathto_seg_wes] is np.nan or wes[pathto_cnvpng_wes] is np.nan: # skip samples that don't have this information
              continue
          if wes[pathto_seg_wes] is not np.nan: # might be smarter to just skip samples without this info, taking note that we're skipping?
              os.system('gsutil cp ' + wes[pathto_seg_wes] + ' ' + tempdir + 'wes_copy_number.tsv')
          if wes[pathto_cnvpng_wes] is not np.nan:
              os.system('gsutil cp ' + wes[pathto_cnvpng_wes] + ' ' + tempdir + external_id + '_wes_copy_number_map.pdf')
          os.system('sips -s format png ' + tempdir + external_id + '_wes_copy_number_map.pdf --out ' + tempdir + external_id + '_wes_copy_number_map.png')
          # annotate image with sample information
          imagedir = tempdir + external_id + '_wes_copy_number_map.png' # will this ever not exist? might want to check if filepath exists
          text = wes["external_id_capture"].replace("_", " ") + " (WES)"
          Helper.addTextToImage(imagedir, text, outputpath=imagedir)
          images.append(tempdir + external_id + '_wes_copy_number_map.png')
          os.system('gsutil cp ' + imagedir + ' ' + outputloc + external_id + '_wes_copy_number_map.png')

          if is_from_pairs:
            snvs = pair_wes[pair_wes["case_sample"] == k]
            if pathto_snv_wes in snvs:
              snvs = snvs[pathto_snv_wes]
            else:
              continue
            for snv in snvs:
              if snv is not np.nan:
                os.system('gsutil cp ' + snv + ' ' + tempdir + 'wes_mutations.tsv')
                # add first matched normal to list if exists (aka add the normal's sample_id)
                try:
                  samples_for_CN_heat_WES += [next(id for id in pair_subset["control_sample"] if id not in ["NA",np.nan])]
                except:
                  continue
                break
          else:
            os.system('gsutil cp ' + wes[pathto_snv_wes] + ' ' + tempdir + 'wes_mutations.tsv')
          # add culture condition and external_id columns to the tables
          wes_cn = pd.read_csv(tempdir + 'wes_copy_number.tsv', sep='\t')
          wes_cn['external_id'] = external_id
          wes_cnfile = cnfile.append(wes_cn)
          wes_mut = pd.read_csv(tempdir + 'wes_mutations.tsv', sep='\t')
          wes_mut['external_id'] = external_id
          wes_mutfile = mutfile.append(wes_mut)
      if found:
        wes_cnfile.to_csv(tempdir + 'wes_cn.tsv', sep='\t')
        wes_mutfile.to_csv(tempdir + 'wes_mut.tsv', sep='\t')
        # wes_ext_ids.to_csv(tempdir + 'wes_all_external_ids.tsv', sep='\t')
        os.system('gsutil cp ' + tempdir + 'wes_cn.tsv ' + outputloc + 'wes_copy_number.tsv')
        os.system('gsutil cp ' + tempdir + 'wes_mut.tsv ' + outputloc + 'wes_mutations.tsv')
        # os.system('gsutil cp ' + tempdir + 'all_external_ids.tsv ' + outputloc + 'wes_all_external_ids.tsv')

      # get CN plot for matched normal if a matched normal exists for a given participant (assuming only 1 distinct matched normal exists)
      matched_normal_ids = pair_wes[pair_wes["case_sample"] == k]["control_sample"]
      for normal_id in matched_normal_ids:
        if normal_id is not np.nan:
          # get the image
          normal_sample = sample_wes[sample_wes.index == normal_id] # indexed by sample id
          if normal_sample[pathto_cnvpng_wes][normal_id] is not np.nan:
            os.system('gsutil cp ' + normal_sample[pathto_cnvpng_wes][normal_id] + ' ' + tempdir + normal_id + '_wes_copy_number_map.png')

            # annotate matched normal CN image
            imagedir = tempdir + normal_id + '_wes_copy_number_map.png'
            if os.path.exists(imagedir):
              text = normal_sample["external_id_capture"][normal_id].replace("_", " ") + " (WES, matched normal)"
              Helper.addTextToImage(imagedir, text, outputpath=imagedir)
              images.append(tempdir + normal_id + '_wes_copy_number_map.png')
              os.system('gsutil cp ' + imagedir + ' ' + outputloc + normal_id + '_wes_copy_number_map.png')

            break
      # # annotate matched normal CN image
      # imagedir = tempdir + normal_id + '_wes_copy_number_map.png'
      # if os.path.exists(imagedir):
      #   text = normal_sample["external_id_capture"][normal_id].replace("_", " ") + " (WES, matched normal)"
      #   Helper.addTextToImage(imagedir, text, outputpath=imagedir)
      #   images.append(tempdir + normal_id + '_wes_copy_number_map.png')
      #   os.system('gsutil cp ' + imagedir + ' ' + outputloc + normal_id + '_wes_copy_number_map.png')

    # add the final list of all samples/external_ids for the given participant
    ext_ids = pd.DataFrame(ext_ids, columns=["external_id"])
    ext_ids.to_csv(tempdir + 'all_external_ids.tsv', sep='\t', index=False)
    os.system('gsutil cp ' + tempdir + 'all_external_ids.tsv ' + outputloc + 'all_external_ids.tsv')

    # # upload all images stored in images to the google bucket
    # for img_path in images:
    #   os.system('gsutil cp ' img_path + ' ' + outputloc + 'something_here.tsv')


    if not found:
      print("we did not find any data for " + val)
    else:
      Helper.mergeImages(images, tempdir + 'merged.png')
      os.system('gsutil cp ' + tempdir + 'merged.png ' + outputloc + 'merged_copy_number_map.png')