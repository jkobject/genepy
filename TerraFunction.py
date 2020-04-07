# Jeremie Kalfon
# jkobject@gmail.com
# for BroadInstitute CDS
# in 2019

import time
import pandas as pd
import dalmatian as dm
import numpy as np
import os
import re
import signal
from JKBio import Helper as h
from JKBio import GCPFunction as gcp
import ipdb


def createManySubmissions(workspace, workflow, references, entity=None, expression=None, use_callcache=True):
  """
  wrapper to create many submissions for a workflow

  Args:
  ----
    workspace: str namespace/workspace from url typically
      namespace (str): project to which workspace belongs
      workspace (str): Workspace name
    references: list(str) a list of name of the row in this entity
    entity: str terra csv type (sample_id...)
    expresson: str to use if want to compute on the direct value of the entity or on values of values
                e.g. this.samples
    use_callcache: Bool to false if want to recompute everything even if same input

  Returns:
  ------
    submission_ids list(str) the submission ids
  """
  wm = dm.WorkspaceManager(workspace)
  submission_ids = []
  for ref in references:
    submission_ids += [wm.create_submission(workflow, ref, entity, expression, use_callcache)]
  return submission_ids


def waitForSubmission(workspace, submissions, raise_errors=True):
  """
  wrapper to create many submissions for a workflow

  Args:
  -----
    workspace: str namespace/workspace from url typically
      namespace (str): project to which workspace belongs
      workspace (str): Workspace name
    submissions: list[str] of submission ids
    raise_errors: to true if errors should stop your code

  Returns:
  -------

  Raises:
  ------
  """
  failed_submission = []
  timing = 0
  wm = dm.WorkspaceManager(workspace)
  assert submissions is not None
  if type(submissions) is type(""):
    submissions = [submissions]
  for scount, submission_id in enumerate(submissions):
    finished = False
    while not finished:
      done = 0
      failed = 0
      finished = True
      submission = wm.get_submission(submission_id)["workflows"]
      for wcount, i in enumerate(submission):
        if i['status'] not in {'Done', 'Aborted', 'Failed', 'Succeeded'}:
          finished = False
        elif i["status"] in {'Failed', 'Aborted'}:
          failed += 1
          if i["workflowEntity"]["entityName"] not in failed_submission:
            print(i["workflowEntity"]["entityName"])
            failed_submission.append(i["workflowEntity"]["entityName"])
        elif i["status"] in {'Done', 'Succeeded'}:
          done += 1
      if not finished:
        time.sleep(40)
        print("status is: Done for " + str(done) + " jobs in submission " + str(scount) + ". " + str(timing) + ",5 mn elapsed.", end="\r")
        timing += 1
        time.sleep(20)
        print("status is: Failed for " + str(failed) + " jobs in submission " + str(scount) + ". " + str(timing) + " mn elapsed.", end="\r")
      else:
        print(str(done / (done + failed)) + " of jobs Succeeded in submission " + str(scount) + ".")
  if len(failed_submission) > 0 and raise_errors:
    raise RuntimeError(str(len(failed_submission)) + " failed submission")
  return failed_submission
  # print and return well formated data


def uploadFromFolder(gcpfolder, prefix, workspace, sep='_', updating=False,
                     fformat="fastq12", newsamples=None, samplesetname=None, source='U'):
  """
  upload samples (virtually: only creates tsv file) from a google bucket to a terra workspace

  A very practical function when you have uploaded a folder of samples (RNA/WES/...) in google storage
  with some naming convention, and want to have them all well listed in Terra for futher processing
  it can create a sample set.
  for a set of files: gs://bucket/path/to/files


  Args:
  -----
    gcpfolder: a gs folder path
    prefix: str the folder path
    workspace: str namespace/workspace from url typically
      namespace (str): project to which workspace belongs
      workspace (str): Workspace name
    sep: str the separator (only takes the first part of the name before the sep character)
    fformat bambai, fastq12, fastqR1R2 given the set of files in the folder (they need to be in this naming format)
            e.g. name.bam name.bai / name1.fastq name2.fastq / name_R1.fastq name_R2.fastq
    newsamples: DONT USE
    samplesetname: str all uploaded samples should be part of a sampleset with name..
  """
  wm = dm.WorkspaceManager(workspace)
  print('please be sure you gave access to your terra email account access to this bucket')
  if samplesetname is None:
    samplesetname = 'from:' + gcpfolder + prefix
  files = gcp.list_blobs_with_prefix(gcpfolder, prefix, '/')
  if fformat == "bambai":
    if newsamples is None:
      data = {'sample_id': [], 'bam': [], 'bai': []}
      for file in files:
        if val.split('.')[-1] in ["bam", "bai"]:
          name = file.split('/')[-1].split('.')[0].split(sep)[0][:-2]
          if name in data['sample_id']:
            pos = data['sample_id'].index(name)
            if file[-4:] == ".bam":
              data['bam'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            elif file[-4:] == ".bai":
              data['bai'].insert(pos, 'gs://' + gcpfolder + '/' + file)
          else:
            data['sample_id'].append(name)
            if file[-4:] == ".bam":
              data['bam'].append('gs://' + gcpfolder + '/' + file)
            elif file[-4:] == ".bai":
              data['bai'].append('gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
        else:
          print("unrecognized file type : " + file)
      df = pd.DataFrame(data)
      df = df.set_index("sample_id")
      df["participant"] = pd.Series(data['sample_id'], index=data['sample_id'])
      wm.upload_samples(df)
      wm.update_sample_set(samplesetname, df.index.values.tolist())
    else:
      # TODO: check if each column exists and can be added, else don't add it
      for i, val in enumerate(newsample["file_path"]):
        if val.split('/')[-1].split('.')[1] != "WholeGenome" or val.split('/')[-2] != "bam":
          newsample = newsample.drop(i)
        elif val.split('/')[1] != 'gs:':
          newsample["file_path"][i] = gcpfolder + newsample["file_path"][i].split('/')[-1]
      newsample = newsample.reset_index(drop=True)
      newsample = newsample.rename(index=str, columns={"sample_name": "sample_id", "subject_name": "participant_id", "file_path": "WGS_bam"})
      currfile = ""
      bai = [''] * int(newsample.shape[0])
      # creating an array of bai and adding it to their coresponding bams
      for i in newsample.index:
        currfile = newsample["WGS_bam"][i]
        if currfile.split('/')[-1].split('.')[-1] == "bai":
          bai[int(newsample[newsample["WGS_bam"] == currfile[:-4]].index.values[0])] = currfile
      newsample["WGS_bam_index"] = pd.Series(bai, index=newsample.index)
      # removing original bai rows
      for i in newsample.index:
        currfile = newsample["WGS_bam"][i]
        if currfile.split('/')[-1].split('.')[-1] == "bai":
          newsample = newsample.drop(i)
      newsample = newsample.reset_index(drop=True)
      newsample["sample_set"] = pd.Series([samplesetname] * int(newsample.shape[0]), index=newsample.index)
      newsample.set_index("sample_id", inplace=True, drop=True)
      newsample = newsample[newsample.columns.tolist()[1:] + [newsample.columns.tolist()[0]]]
      newsample = newsample.loc[~newsample.index.duplicated(keep='first')]
      newsample.to_csv("temp/samples.bambai.tsv", sep="\t")
      wm.upload_samples(newsample)
      wm.update_sample_set(samplesetname, newsample.index)
  if fformat in {"fastq12", "fastqR1R2"}:
    data = {'sample_id': [], 'fastq1': [], 'fastq2': []}
    # print and return well formated data
    print(files)
    for file in files:
      if file[-9:] == ".fastq.gz" or file[-6:] == ".fq.gz":
        name = file.split('/')[-1].split('.')[0].split(sep)[0]
        if name in data['sample_id']:
          pos = data['sample_id'].index(name)
          if fformat == "fastqR1R2":
            if "R1" in file:
              data['fastq1'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            elif "R2" in file:
              data['fastq2'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
          else:
            if file.split('.')[-3][-1] == '1':
              data['fastq1'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            elif file.split('.')[-3][-1] == '2':
              data['fastq2'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq 1/2 error", file)
        else:
          data['sample_id'].append(name)
          if fformat == "fastqR1R2":
            if "R1" in file:
              data['fastq1'].append('gs://' + gcpfolder + '/' + file)
            elif "R2" in file:
              data['fastq2'].append('gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
          else:
            if file.split('.')[-3][-1] == '1':
              data['fastq1'].append('gs://' + gcpfolder + '/' + file)
            elif file.split('.')[-3][-1] == '2':
              data['fastq2'].append('gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
      else:
        print("unrecognized file type : " + file)
    df = pd.DataFrame(data)
    print(df)
    df["Source"] = source
    df["participant"] = data['sample_id']
    df = df.set_index("sample_id")
    wm.upload_samples(df)
    wm.update_sample_set(samplesetname, df.index.values.tolist())


def updateAllSampleSet(workspace, newsample_setname, Allsample_setname='All_samples'):
  """
  update the previous All Sample sample_set with the new samples that have been added.

  It is especially useful for the aggregate task. Can more generally merge two samplesets together

  Args:
  ----
    workspace: str namespace/workspace from url typically
      namespace (str): project to which workspace belongs
      workspace (str): Workspace name
    newsample_setname: str name of sampleset to add to All_samples
  """
  prevsamples = list(dm.WorkspaceManager(workspace).get_sample_sets().loc[Allsample_setname]['samples'])
  newsamples = list(dm.WorkspaceManager(workspace).get_sample_sets().loc[newsample_setname]['samples'])
  prevsamples.extend(newsamples)
  dm.WorkspaceManager(workspace).update_sample_set(Allsample_setname, list(set(prevsamples)))


def addToSampleSet(workspace, samplesetid, samples):
  """

  will create new if doesn't already exist, else adds to existing

  Args:
  ----
  """
  try:
    prevsamples = dm.WorkspaceManager(workspace).get_sample_sets()['samples'][samplesetid]
    samples.extend(prevsamples)
  except KeyError:
    print('The sample set ' + str(samplesetid) + ' did not exist in the workspace. Will be created now...')
  dm.WorkspaceManager(workspace).update_sample_set(samplesetid, list(set(samples)))


def addToPairSet(workspace, pairsetid, pairs):
  """
  Similar to above but for a pairset
  """

  try:
    prevpairs = dm.WorkspaceManager(workspace).get_pair_sets().loc[[pairsetid]].pairs[0]
    pairs.extend(prevpairs)
  except KeyError:
    print('The pair set ' + str(pairsetid) + ' did not exist in the workspace. Will be created now...')
  dm.WorkspaceManager(workspace).update_pair_set(pairsetid, list(set(pairs)))

# Gwen's old version - caught some niche conditions made by get_pair_sets()
# that I think may raise errors in the current version.
# yep. the niche condition has to do with the fact that a list isn't hashable
# def addToPairSet(wm, pairsetid, pairs):
#   pairsets = wm.get_pair_sets()
#   prevpairs = pairsets.loc[[pairsetid]].pairs.tolist() # is this always a list of list? I think so.
#   print(type(prevpairs[0]))
#   if isinstance(prevpairs[0], str) :
#     pairs.extend(prevpairs)
#   elif isinstance(prevpairs[0], list):
#     pairs.extend(prevpairs[0])
#   wm.update_pair_set(pairsetid, list(set(pairs)))


def saveOmicsOutput(workspace, pathto_cnvpng='segmented_copy_ratio_img',
                    pathto_stats='sample_statistics',
                    specific_cohorts=[],
                    speicifc_celllines=[],
                    is_from_pairs=True,
                    pathto_snv='filtered_variants',
                    pathto_seg='cnv_calls',
                    datadir='gs://cclf_results/targeted/kim_sept/',
                    specific_samples=[]):
  """

  """
  if specific_cohorts:
    samples = dm.WorkspaceManager(workspace).get_samples()
    samples = samples[samples.index.isin(specificlist)]
  if is_from_pairs:
    pairs = dm.WorkspaceManager(workspace).get_pairs()
    pairs = pairs[pairs['case_sample'].isin(specificlist)]
  for i, val in samples.iterrows():
    os.system('gsutil cp ' + val[pathto_seg] + ' ' + datadir + i + '/')
    os.system('gsutil cp ' + val[pathto_cnvpng] + ' ' + datadir + i + '/')
    os.system('gsutil cp ' + val[pathto_stats] + ' ' + datadir + i + '/')
    if is_from_pairs:
      snvs = pairs[pairs["case_sample"] == i][pathto_snv]
      for snv in snvs:
        if snv is not np.nan:
          os.system('gsutil cp ' + snv + ' ' + datadir + i + '/')
          break
    else:
      os.system('gsutil cp ' + val[pathto_snv] + ' ' + datadir + i + '/')


def changeGSlocation(workspacefrom, workspaceto=None, prevgslist=[], newgs='', index_func=None,
                     flag_non_matching=False, onlycol=[], entity='', droplists=True, keeppath=True, dry_run=False):
  """

  """
  flaglist = []
  data = {}
  wmfrom = dm.WorkspaceManager(workspacefrom)
  if not keeppath:
    flag_non_matching = True
  if entity in ['', 'participants']:
    try:
      a = wmfrom.get_participants()
      data.update({'participants': a})
    except:
      print('no participants')
  if entity in ['', 'samples']:
    try:
      a = wmfrom.get_samples()
      data.update({'samples': a})
    except:
      print('no samples')
  if entity in ['', 'pair_sets']:
    try:
      a = wmfrom.get_pair_sets()
      data.update({'pair_sets': a})
    except:
      print('no pair_sets')
  if entity in ['', 'pairs']:
    try:
      a = wmfrom.get_pairs()
      data.update({'pairs': a})
    except:
      print('no pairs')
  if entity in ['', 'sample_sets']:
    try:
      a = wmfrom.get_sample_sets()
      data.update({'sample_sets': a})
    except:
      print('no sample_sets')
    # currently works only for sample, sample
  for i, entity in data.items():
    if onlycol:
      try:
        entity = entity[onlycol]
      except:
        print("entity " + str(i) + " does not contain one of the columns")
        continue
    todrop = set()
    for j, val in entity.iterrows():
      # print(j)
      for k, prev in enumerate(val):
        if type(prev) is str:
          new = prev
          if newgs not in new:
            for prevgs in prevgslist:
              new = new.replace(prevgs, newgs)
            if flag_non_matching:
              if 'gs://' == prev[:5]:
                if new == prev:
                  flaglist.append(prev)
                elif not keeppath:
                  new = newgs + new.split('/')[-1]
          val[k] = new
        if type(prev) is list:
          if droplists:
            todrop.add(k)
            continue
          ind = []
          for prevname in prev:
            newname = prevname
            if newgs not in newname:
              for prevgs in prevgslist:
                newname = newname.replace(prevgs, newgs)
              if flag_non_matching:
                if 'gs://' == prevname[:5]:
                  if newname == prevname:
                    flaglist.append(prevname)
                  elif not keeppath:
                    new = newgs + new.split('/')[-1]
            ind.append(newname)
          val[k] = ind
        entity.loc[j] = val
    if onlycol:
      data[i][onlycol] = entity
    for drop in todrop:
      data[i] = data[i].drop(drop, 1)
  if workspaceto is None:
    wmto = wmfrom
  else:
    wmto = dm.WorkspaceManager(workspaceto)
  for key in data.keys():
    for k in data[key].columns:
      data[key][k] = data[key][k].astype(str)
  if not dry_run:
    ipdb.set_trace()
    if "participants" in data:
      wmto.upload_entities('participant', data['participants'])
    if "samples" in data:
      wmto.upload_samples(data['samples'])
    if "pairs" in data:
      wmto.upload_pairs(data['pairs'])
    if "pair_set" in data:
      pairset = data['pair_set'].drop('pairs', 1)
      wmto.upload_entities('pair_set', pairset)
    if "sample_sets" in data:
      sampleset = data['sample_sets'].drop('samples', 1)
      wmto.upload_entities('sample_set', sampleset)
  return flaglist


def renametsvs(workspace, wmto=None, index_func=None):
  """
  ################## WIP ############
  only works for one use case
  """
  data = {}
  wmfrom = dm.WorkspaceManager(workspace)
  try:
    a = wmfrom.get_participants()
    data.update({'participants': a})
  except:
    print('no participants')
  try:
    a = wmfrom.get_samples()
    data.update({'samples': a})
  except:
    print('no samples')
  try:
    a = wmfrom.get_pair_sets()
    data.update({'pair_sets': a})
  except:
    print('no pair_sets')
  try:
    a = wmfrom.get_pairs()
    data.update({'pairs': a})
  except:
    print('no pairs')
  try:
    a = wmfrom.get_sample_sets()
    data.update({'sample_sets': a})
  except:
    print('no sample_sets')
  # currently works only for sample, sample
  for k, entity in data.items():
    ind = []
    for i in entity.index:
      pos = val.find('-SM')
      if pos != -1:
        val = val[pos + 1:]
        pos = val.find('-SM')
        if pos != -1:
          val = val[:9] + val[pos + 1:]
      ind.append(val)
    entity.index = ind
    # for all columns of the tsv
    for k, val in entity.iterrows():
      for i, v in enumerate(val):
        if type(v) is list or type(v) is str:
          ind = []
          for j in v:
            pos = j.find('-SM')
            if pos != -1:
              j = j[pos + 1:]
              pos = j.find('-SM')
              if pos != -1:
                j = j[:9] + j[pos + 1:]
            ind.append(j)
          val[i] = ind
        entity.loc[k] = val
    if wmto is None:
      wmto = wmfrom
    if "participants" in data:
      wmto.upload_participants(data['participants'].index.tolist())
    if "samples" in data:
      wmto.upload_samples(data['samples'])
    if "pairs" in data:
      wmto.upload_entities('pair', data['pairs'])
    if "pair_set" in data:
      pairset = data['pair_set'].drop('pairs', 1)
      wmto.upload_entities('pair_set', pairset)
      for i, val in data['pair_set'].iterrows():
        wmto.update_pair_set(i, val.pairs)
    if "sample_set" in data:
      sampleset = data['sample_set'].drop('samples', 1)
      wmto.upload_entities('sample_set', sampleset)
      for i, val in data['sample_set'].iterrows():
        wmto.update_sample_set(i, val.samples)


def findBackErasedDuplicaBamteFromTerraBucket(workspace, gsfolder, bamcol="WES_bam", baicol="WES_bai"):
  """
  If you have erased bam files in gcp with bai files still present and the bam files are stored elsewhere
  and their location is in a terra workspace.

  Will find them back by matching bai sizes and copy them back to their original locations

  Args:
  ----
    workspace: str namespace/workspace from url typically
      namespace (str): project to which workspace belongs
      workspace (str): Workspace name
    gsfolder: str the gsfolder where the bam files are
  """
  # get ls of all files folder
  samples = os.popen('gsutil -m ls -al ' + gsfolder + '**.bai').read().split('\n')
  # compute size filepath

  sizes = {'gs://' + val.split('gs://')[1].split('#')[0]: int(val.split("2019-")[0]) for val in samples[:-2]}
  names = {}
  for k, val in sizes.items():
    if val in names:
      names[val].append(k)
    else:
      names[val] = [k]
  # get all bai in tsv
  samp = dm.WorkspaceManager(workspace).get_samples()
  for k, val in samp.iterrows():
    if val[bamcol] != 'NA' and val[baicol] != 'NA':
      # if bai has duplicate size
      code = os.system('gsutil ls ' + val[bamcol])
      if code == 256:
        if val[bamcol] is None:
          print('we dont have bam value for ' + str(k))
          continue
        else:
          print('no match values for ' + str(val[bamcol]))

        for va in names[sizes[val[baicol]]]:
          # for all duplicate size
          # if ls bam of bai duplicate size work
          # mv bam to bampath in gsfolder
          if '.bam' in va:
            if os.system('gsutil ls ' + va.split('.bam.bai')[0] + ".bam") == 0:
              print('gsutil mv ' + va.split('.bam.bai')[0] + ".bam " + val[bamcol])
              os.system('gsutil mv ' + va.split('.bam.bai')[0] + ".bam " + val[bamcol])
              break
          elif os.system('gsutil ls ' + va.split('.bai')[0] + ".bam") == 0:
            print('gsutil mv ' + va.split('.bai')[0] + ".bam " + val[bamcol])
            os.system('gsutil mv ' + va.split('.bai')[0] + ".bam " + val[bamcol])
            break
      elif code == signal.SIGINT:
        print('Awakened')
        break
    else:
      print("no data for " + str(k))


def shareTerraBams(users, workspace, samples, bamcols=["WES_bam", "WES_bai"]):
  """
  will share some files from gcp with a set of users using terra as metadata repo.

  only works with files that are listed on a terra workspace tsv but actually
  point to a regular google bucket and not a terra bucket.

  Args:
  ----
    users: list[str] of users' google accounts
    workspace: str namespace/workspace from url typically
      namespace (str): project to which workspace belongs
      workspace (str): Workspace name
    samples list[str] of samples_id for which you want to share data
    bamcols: list[str] list of column names
  """
  if type(users) is str:
    users = [users]
  wm = dm.WorkspaceManager(workspace)
  togiveaccess = np.ravel(wm.get_samples()[bamcols].loc[samples].values)
  for user in users:
    files = ''
    for i in togiveaccess:
      files += ' ' + i
    code = os.system("gsutil -m acl ch -ru " + user + ":R" + files)
    if code == signal.SIGINT:
      print('Awakened')
      break
  print('the files are stored here:\n\n')
  print(togiveaccess)
  print('\n\njust install and use gsutil to copy them')
  print('https://cloud.google.com/storage/docs/gsutil_install')
  print('https://cloud.google.com/storage/docs/gsutil/commands/cp')
  return togiveaccess


def saveConfigs(workspace, filepath):
  """
  will save everything about a workspace into a csv and json file

  Args:
  -----
    workspace: str namespace/workspace from url typically
      namespace (str): project to which workspace belongs
      workspace (str): Workspace name
    filepath to save files
  """
  wm = dm.WorkspaceManager(workspace)
  h.createFoldersFor(filepath)

  conf = wm.get_configs()
  conf.to_csv(filepath + '.csv')
  params = {}
  params['GENERAL'] = wm.get_workspace_metadata()
  for k, val in conf.iterrows():
    params[k] = wm.get_config(val['name'])
  h.dictToFile(params, filepath + '.json')
