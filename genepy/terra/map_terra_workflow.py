# 1. make sure dalmatian is installed
# pip install firecloud-dalmatian
# 2. make sure graphviz is installed (this script will run "dot" to generate images)
# brew install graphviz

import subprocess
import requests
import pandas as pd
import dalmatian
import subprocess


def resolve_dot_path(root_entity_type, name):
  name = name.replace(" ", "")
  parts = name.split(".")
  assert parts[0] == "this"
  del parts[0]
  cur = root_entity_type
  while len(parts) > 1:
    next_part = parts[0]
    if cur == "sample_set" and next_part == "samples":
      cur = "sample"
    elif cur == "pair" and next_part == "case_sample":
      cur = "sample"
    elif cur == "pair" and next_part == "control_sample":
      cur = "sample"
    elif cur == "sample" and next_part == "participant":
      cur = "participant"
    else:
      raise Exception(f"Unknown case: {cur} -> {next_part} (root_entity_type={root_entity_type} name={name})")
    del parts[0]
  return cur+"."+parts[0]

def extract_config_summary(workspace_name, workflows=None):
  wm = dalmatian.WorkspaceManager(workspace_name)
  configs = wm.get_configs()

  config_summaries = []
  for rec in configs.to_records():
    cfgname = rec['namespace']+"/"+rec['name']
    if workflows is not None:
      if cfgname not in workflows:
        continue
    config = wm.get_config(cfgname)
    config['inputs'] = {k:v.strip() for k,v in config['inputs'].items()}
    config['outputs'] = {k:v.strip() for k,v in config['outputs'].items()}
    inputs=[resolve_dot_path(config['rootEntityType'], x) for x in config['inputs'].values() if x.startswith("this.")]
    outputs=[resolve_dot_path(config['rootEntityType'], x) for x in config['outputs'].values() if x.startswith("this.")]
    config_summaries.append(dict(inputs=inputs, outputs=outputs, entity_type=rec['rootEntityType'], name=cfgname))
  return config_summaries


def write_dependency_graph_image(filename, config_summaries):
  with open("/tmp/sample.dot", "wt") as fd:
    node_names = {}
    def nn(name, is_var):
      if name in node_names:
        return node_names[name]
      node_name = "n{}".format(len(node_names))
      node_names[name] = node_name
      fd.write("{} [label=\"{}\" {}];\n".format(node_name, name, {True: "shape=oval", False: "shape=box fillcolor=yellow style=filled"}[is_var]))
      return node_names[name]

    fd.write("digraph { rankdir=LR;\n")
    for config in config_summaries:
      for name in config['inputs']:
        fd.write("{} -> {};\n".format(nn(name, True), nn(config['name'], False)))

      for name in config['outputs']:
        fd.write("{} -> {};\n".format(nn(config['name'], False), nn(name, True)))
    fd.write("}\n")


  subprocess.check_call(["dot", "/tmp/sample.dot", "-Tpng", "-o", filename])

def write_config_summary_table(filename, config_summaries):
  from collections import defaultdict
  variables = defaultdict(lambda: dict(used_by=[], produced_by=[]))
  for config in config_summaries:
    for name in config['inputs']:
      variables[name]['used_by'].append(config['name'])

    for name in config['outputs']:
      variables[name]['produced_by'].append(config['name'])

  var_col = []
  used_by = []
  produced_by = []
  inputs = []
  for k, v in variables.items():
    if len(v["produced_by"]) == 0:
      inputs.append(k)
    var_col.append(k)
    used_by.append(", ".join(v["used_by"]))
    produced_by.append(", ".join(v["produced_by"]))
  df = pd.DataFrame(dict(variable=var_col, used_by=used_by, produced_by=produced_by))
  df.to_csv(filename)

def map_workspace_diagram(workspace_name, output_path='terra-workflows', workflows=None):
  """
  -- adapted from scripts written by Philip Montgomery --
  this function creates a graph of the workflows within a Terra workspace

  inputs:
      workspace_name (str): name of the workspace
      output_path (str): path where outputs will be saved
      workflows (list[str]): list of workflows to consider in the graph. If None is provided (default),
          all the workflows in the workspace will be used.

  example code:
      workflows = ['stewart/pipette_wgs_SV', 'stewart/manta', 'stewart/SvABA_xtramem',
              'stewart/svaba_snowmanvcf2dRangerForBP', 'stewart/mantavcf2dRangerForBP',
              'stewart/extract_dRanger_intermediates','stewart/pcawg_snowmanvcf2dRangerForBP',
              'stewart/SV_cluster_forBP', 'stewart/breakpointer', 'stewart/Breakpointer_fix_sample',
              'stewart/REBC_SV_consensus_filter_v3']

      workspace_name = 'broad-firecloud-ccle/REBC_methods_only-tmp'
      print(workspace_name)
      config = mtw.map_workspace_diagram(workspace_name, workflows=workflows)
  """
  configs = extract_config_summary(workspace_name, workflows=workflows)
  write_dependency_graph_image(output_path+'/'+workspace_name.replace("/", " ")+".png", configs)
  write_config_summary_table(output_path+'/'+workspace_name.replace("/", " ")+".csv", configs)
  return configs
