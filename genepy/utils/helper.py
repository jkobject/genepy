# Jeremie Kalfon
# for BroadInsitute
# in 2019

from __future__ import print_function
import json
import os
import string
import subprocess
import pandas as pd
import numpy as np
import itertools
import random
import xmltodict
from biomart import BiomartServer
import io

chromsize = {
    "chr1": 248956422,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr2": 242193529,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chrX": 156040895,
    "chrY": 57227415,
}

rename_mut = {
    "contig": "chr",
    "position": "pos",
    "Reference_Allele": "ref",
    "ref_allele": "ref",
    "alt_allele": "alt",
    "Chromosome": "chr",
    "End_postition": "end",
    "Start_position": "pos",
    "Tumor_Seq_Allele1": "alt",
}


def fileToList(filename, strconv=lambda x: x):
    """
    loads an input file with a\\n b\\n.. into a list [a,b,..]
    """
    with open(filename) as f:
        return [strconv(val[:-1]) for val in f.readlines()]


def listToFile(l, filename, strconv=lambda x: str(x)):
    """
    loads a list with [a,b,..] into an input file a\\n b\\n..
    """
    with open(filename, "w") as f:
        for item in l:
            f.write("%s\n" % strconv(item))


def dictToFile(d, filename):
    """
    turn a dict into a json file
    """
    with open(filename, "w") as json_file:
        json.dump(d, json_file)


def fileToDict(filename):
    """
    loads a json file into a python dict
    """
    with open(filename) as f:
        data = json.load(f)
    return data


def batchMove(l, pattern=["*.", ".*"], folder="", add=""):
    """
    moves a set of files l into a folder:

    Args:
      l: file list
      pattern: if files are a set of patterns to match
      folder: folder to move file into
      add: some additional mv parameters
    """
    for val in l:
        cmd = "mv "
        if add:
            cmd += add + " "
        if "*." in pattern:
            cmd += "*"
        cmd += val
        if ".*" in pattern:
            cmd += "*"
        cmd += " " + folder
        res = os.system(cmd)
        if res != 0:
            raise Exception("Leave command pressed or command failed")


def batchRename(dt, folder="", sudo=False, doAll=False, add="", dryrun=False):
    """
    Given a dict renames corresponding files in a folder

    Args:
      dt (dict): dict(currentName:newName) renaming dictionnary
      folder (str): folder to look into
      add: some additional mv parameters
    """
    cmd = "ls -R " + folder if doAll else "ls " + folder
    files = os.popen(cmd).read().split("\n")
    if doAll:
        prep = ""
        f = []
        for val in files:
            if len(val) == 0:
                prep = ""
                continue
            if val[0] == "." and len(val) > 3:
                prep = val[:-1]
                continue
            if "." in val:
                f.append(prep + "/" + val)
        files = f
    for k, val in dt.items():
        for f in files:
            if k in f:
                cmd = "sudo mv " if sudo else "mv "
                if add:
                    cmd += add + " "
                if not doAll:
                    cmd += folder
                cmd += f
                cmd += " "
                if not doAll:
                    cmd += folder
                cmd += f.replace(k, val)
                if dryrun:
                    print(cmd)
                else:
                    res = os.system(cmd)
                    if res != 0:
                        raise Exception("Leave command pressed or command failed")


def grouped(iterable, n):
    """
    iterate over element of list 2 at a time python

    s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ...
    """
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk


def overlap(interval1, interval2):
    """computed overlap

    Given [0, 4] and [1, 10] returns [1, 4]
    Given [0, 4] and [8, 10] returns False
    """
    if interval2[0] <= interval1[0] <= interval2[1]:
        start = interval1[0]
    elif interval1[0] <= interval2[0] <= interval1[1]:
        start = interval2[0]
    else:
        return False

    if interval2[0] <= interval1[1] <= interval2[1]:
        end = interval1[1]
    elif interval1[0] <= interval2[1] <= interval1[1]:
        end = interval2[1]
    else:
        return False

    return (start, end)


def union(interval1, interval2):
    """
    Given [0, 4] and [1, 10] returns [0, 10]
    Given [0, 4] and [8, 10] returns False
    """
    if interval1[0] <= interval2[0] <= interval1[1]:
        start = interval1[0]
        end = interval1[1] if interval2[1] <= interval1[1] else interval2[1]
    elif interval1[0] <= interval2[1] <= interval1[1]:
        start = interval2[0] if interval2[0] <= interval1[0] else interval1[0]
        end = interval1[1]
    else:
        return False
    return (start, end)


def nans(df):
    return df[df.isnull().any(axis=1)]


def createFoldersFor(filepath):
    """
    will recursively create folders if needed until having all the folders required to save the file in this filepath
    """
    prevval = ""
    for val in os.path.expanduser(filepath).split("/")[:-1]:
        prevval += val + "/"
        if not os.path.exists(prevval):
            os.mkdir(prevval)


def randomString(stringLength=6, stype="all", withdigits=True):
    """
    Generate a random string of letters and digits

    Args:
      stringLength (int, optional): the amount of char. Defaults to 6.
      stype (str, optional): one of lowercase, uppercase, all. Defaults to 'all'.
      withdigits (bool, optional): digits allowed in the string? Defaults to True.

    Returns:
      str: random string
    """
    if stype == "lowercase":
        lettersAndDigits = string.ascii_lowercase
    elif stype == "uppercase":
        lettersAndDigits = string.ascii_uppercase
    else:
        lettersAndDigits = string.ascii_letters
    if withdigits:
        lettersAndDigits += string.digits
    return "".join(random.choice(lettersAndDigits) for i in range(stringLength))


def pdDo(df, op="mean", of="value1", over="value2"):
    """
    apply a function to a panda dataframe WIP
    """
    df = df.sort_values(by=over)
    index = []
    data = df.iloc[0, of]
    ret = []
    prev = df.iloc[0, over]
    j = 0
    for k, val in df.iloc[1:].iterrows():
        if val[over] == prev:
            data.append(val[of])
        else:
            if of == "mean":
                ret[j] = np.mean(data)
            elif of == "sum":
                ret[j] = np.sum(data)
            elif of == "max":
                ret[j] = np.max(data)
            elif of == "min":
                ret[j] = np.min(data)
            index.append(k)
            j += 1
            data = [val[of]]
    return index, ret


def parrun(cmds, cores, add=[]):
    """
    runs a set of commands in parallel using the "&" command

    Args:
      cmds: the list of commands
      cores: number of parallel execution
      add: an additional list(len(cmds)) of command to run in parallel at the end of each parallel run
    """
    count = 0
    exe = ""
    if len(add) != 0 and len(add) != len(cmds):
        raise ValueError("we would want them to be the same size")
    else:
        addexe = ""
    fullres = []
    for i, cmd in enumerate(cmds):
        count += 1
        exe += cmd
        if len(add) != 0:
            addexe += add[i]
        if count < cores and i < len(cmds) - 1:
            exe += " & "
            if len(add) != 0:
                addexe += " & "
        else:
            count = 0
            res = subprocess.run(exe, capture_output=True, shell=True)
            if res.returncode != 0:
                raise ValueError("issue with the command: " + str(res.stderr))
            exe = ""
            if len(add) != 0:
                res = subprocess.run(addexe, capture_output=True, shell=True)
                if res.returncode != 0:
                    raise ValueError("issue with the command: " + str(res.stderr))
                addexe = ""
            fullres.append(res.stdout.decode("utf-8"))
    return fullres


def askif(quest):
    """
    asks a y/n question to the user about something and returns true or false given his answer
    """
    print(quest)
    inp = input()
    if inp in ["yes", "y", "Y", "YES", "oui", "si"]:
        return 1
    elif inp in ["n", "no", "nope", "non", "N"]:
        return 0
    else:
        return askif("you need to answer by yes or no")


def inttodate(i, lim=1965, unknown="U", sep="-", order="asc", startsatyear=0):
    """
    transforms an int representing days into a date

    Args:
      i: the int
      lim: the limited year below which we have a mistake
      unknown: what to return when unknown (date is bellow the limited year)
      sep: the sep between your date (e.g. /, -, ...)
      order: if 'asc', do d,m,y else do y,m,d
      startsatyear: when is the year to start counting for this int

    Returns:
      str: the date or unknown
    """
    a = int(i // 365)
    if a > lim:
        a = str(a + startsatyear)
        r = i % 365
        m = str(int(r // 32)) if int(r // 32) > 0 else str(1)
        r = r % 32
        d = str(int(r)) if int(r) > 0 else str(1)
    else:
        return unknown
    return d + sep + m + sep + a if order == "asc" else a + sep + m + sep + d


def datetoint(dt, split="-", unknown="U", order="des"):
    """
    same as inttodate but in the opposite way;

    starts at 0y,0m,0d

    Args:
      dt: the date string
      split: the splitter in the string (e.g. /,-,...)
      unknown: maybe the some dates are 'U' or 0 and the program will output 0 for unknown instead of crashing
      order: if 'asc', do d,m,y else do y,m,d

    Returns:
      int: the date
    """
    arr = np.array(dt[0].split(split) if dt[0] != unknown else [0, 0, 0]).astype(int)
    if len(dt) > 1:
        for val in dt[1:]:
            arr = np.vstack(
                (
                    arr,
                    np.array(
                        val.split(split)
                        if val != unknown and val.count(split) == 2
                        else [0, 0, 0]
                    ).astype(int),
                )
            )
        arr = arr.T
    res = (
        arr[2] * 365 + arr[1] * 31 + arr[0]
        if order == "asc"
        else arr[0] * 365 + arr[1] * 31 + arr[2]
    )
    return [res] if type(res) is np.int64 else res


prevshowcount = 100


def showcount(i, size):
    """
    pretty print of i/size%, to put in a for loop
    """
    global prevshowcount
    a = 1 + int(100 * (i / size))
    if a != prevshowcount:
        print(str(a) + "%", end="\r")
        prevshowcount = a


def combin(n, k):
    """
    Nombre de combinaisons de n objets pris k a k
    outputs the number of comabination of n object taken k at a time
    """
    if k > n // 2:
        k = n - k
    x = 1
    y = 1
    i = n - k + 1
    while i <= n:
        x = (x * i) // y
        y += 1
        i += 1
    return x


def dups(lst):
    """
    shows the duplicates in a list
    """
    seen = set()
    # adds all elements it doesn't know yet to seen and all other to seen_twice
    seen_twice = set(x for x in lst if x in seen or seen.add(x))
    # turn the set into a list (as requested)
    return list(seen_twice)


def makeCombinations(size, proba):
    """
    produces probability of X event happening at the same time

    pretty usefull for cobinding analysis. wil compute it
    given binomial probabilities of each event occuring and the number of trials

    Args:
      size: int number of trials
      proba: list[float] probabilities of each event occuring
    """
    sums = {i: 0 for i in range(1, size)}
    for i in range(size - 1, 0, -1):
        print(i)
        if sums[i] > 0:
            continue
        print(combin(size + 3, i))
        v = 0
        for j in itertools.combinations(proba, i):
            v += np.prod(j)
        sums[i] = v
    for i in range(size - 1, 0, -1):
        for j in range(i + 1, size):
            icomb = combin(j, i)
            sums[i] -= icomb * sums[j]
    sums[0] = 1 - sum(list(sums.values()))
    return sums


def closest(lst, K):
    """
    returns the index of the value closest to K in a lst
    """
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


def compareDfs(df1, df2):
    """
    compares df1 to df2

    shows col difference, index difference, nans & 0s differences
    """
    nmissmatchCols = set(df1.columns) - set(df2.columns)
    omissmatchCols = set(df2.columns) - set(df1.columns)
    nmissmatchInds = set(df1.index) - set(df2.index)
    omissmatchInds = set(df2.index) - set(df1.index)
    newNAs = df1.isna().sum().sum() - df2.isna().sum().sum()
    new0s = (df1 == 0).sum().sum() - (df2 == 0).sum().sum()
    print("FOUND missmatch Columns NOT IN df2: " + str(nmissmatchCols))
    print("FOUND missmatch Columns NOT IN df1: " + str(omissmatchCols))
    print("FOUND missmatch Index NOT IN df2: " + str(nmissmatchInds))
    print("FOUND missmatch Index NOT IN df1: " + str(omissmatchInds))
    print("FOUND new NAs in df1: " + str(newNAs))
    print("FOUND new 0s in df1: " + str(new0s))
    return nmissmatchCols, omissmatchCols, nmissmatchInds, omissmatchInds, newNAs, new0s


def stringifydict(res):
    """[summary]

    Args:
        res ([type]): [description]

    Returns:
        [type]: [description]
    """
    a = {}
    for k, v in res.items():
        if type(v) is dict:
            a[k] = stringifydict(v)
        else:
            a[str(k)] = v
    return a


def readXMLs(folder=None, file=None, rename=None):
    """[summary]

    Args:
        folder ([type], optional): [description]. Defaults to None.
        file ([type], optional): [description]. Defaults to None.
        rename ([type], optional): [description]. Defaults to None.

    Raises:
        ValueError: [description]

    Returns:
        [type]: [description]
    """
    if file is not None:
        if type(file) is str:
            print("reading 1 file")
            files = [file]
        else:
            print("reading files")
            files = file
    elif folder is not None:
        print("reading from folder")
        files = [i for i in os.listdir(folder) if i.endswith(".xml")]
    else:
        raise ValueError("need folder or file")
    df = pd.DataFrame()
    for file in files:
        res = []
        a = open(file, "r").read()
        xmldict = xmltodict.parse(a)
        data = xmldict["Workbook"]["Worksheet"]["Table"]["Row"]
        for val in data:
            res.append(
                [
                    v["Data"]["#text"] if "#text" in v["Data"] else None
                    for v in val["Cell"]
                ]
            )
        res = pd.DataFrame(data=res[2:], columns=res[0])
        df = df.append(res)
    if rename is not None:
        df = df.rename(columns=rename)
    return df


def makeCellosaurusExport(
    ftp="https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt",
    reload=False,
    dropped=[
        "ACH-002260",
        "ACH-001741",
        "ACH-001189",
    ],
    ancestry={
        "African": [],
        "Native American": [],
        "East Asian, North": [],
        "East Asian, South": [],
        "South Asian": [],
        "European, North": [],
        "European, South": [],
    },
):
    """
    make a df from cellosaurus' human cancer cell line data

    Args:
      ftp (str): the ftp link to the cellosaurus latest txt db export

    Returns:
      pd.DataFrame: indexed by cellosaurus id and with "depmap_id", "id", "disease", "age",
      "sex", "patient_id", "parent_id", "date", "synonyms", "has_issues", "comments"
    """
    print("make sure to have wget installed")
    if reload or not os.path.exists("/tmp/cellosaurus.txt"):
        subprocess.run("wget " + ftp + " -O /tmp/cellosaurus.txt")
    l = open("/tmp/cellosaurus.txt", "r")
    print(
        """

    ---------  ---------------------------     ----------------------
    Line code  Content                         Occurrence in an entry
    ---------  ---------------------------     ----------------------
    ID         Identifier (cell line name)     Once; starts an entry
    AC         Accession (CVCL_xxxx)           Once
    AS         Secondary accession number(s)   Optional; once
    SY         Synonyms                        Optional; once
    DR         Cross-references                Optional; once or more
    RX         References identifiers          Optional: once or more
    WW         Web pages                       Optional; once or more
    CC         Comments                        Optional; once or more
    ST         STR profile data                Optional; once or more
    DI         Diseases                        Optional; once or more
    OX         Species of origin               Once or more
    HI         Hierarchy                       Optional; once or more
    OI         Originate from same individual  Optional; once or more
    SX         Sex of cell                     Optional; once
    AG         Age of donor at sampling        Optional; once
    CA         Category                        Once
    DT         Date (entry history)            Once
    //         Terminator                      Once; ends an entry

  """
    )
    v = []
    CL = False
    Homo = False
    age = "U"
    sex = "U"
    date = "U"
    parent = ""
    ID = ""
    SY = ""
    DI = ""
    atcclink = ""
    dsmzlink = ""
    doublingt = ""
    hasebv = False
    origin = ""
    isMeta = False
    characteristics = ""
    ancestry = ""
    issue = ""
    instability = ""
    transfected = ""
    depmap = set()
    lid = None
    individual = "ID-" + randomString(stringLength=6, stype="all", withdigits=True)
    comments = [""]
    tofind = {}
    cl = {}
    while True:
        # Get next line from file
        line = l.readline()
        # if line is empty
        # end of file is reached
        if not line:
            break
        else:
            line = line[:-1]
        if line[:2] == "AC":
            lid = line[5:]
            if lid in cl:
                lid = None
        elif line[:2] == "//":
            depmap -= set(dropped)
            if ((CL and Homo and ID) or depmap) and (lid is not None):
                if len(depmap) > 1:
                    raise ValueError("multiple depmap mapping to one line", depmap)
                elif len(depmap) == 1:
                    depmap = depmap.pop()
                else:
                    depmap = ""
                v.append(depmap)
                v.append(ID)
                v.append(DI)
                v.append(age)
                v.append(sex)
                v.append(individual)
                v.append(parent)
                v.append(date)
                v.append(SY)
                v.append(" | ".join(comments))
                v.append(atcclink)
                v.append(dsmzlink)
                v.append(doublingt)
                v.append(hasebv)
                v.append(origin)
                v.append(isMeta)
                v.append(characteristics)
                v.append(ancestry)
                v.append(issue)
                v.append(instability)
                v.append(transfected)
                cl.update({lid: v})
            v = []
            CL = False
            Homo = False
            age = ""
            sex = ""
            date = ""
            parent = ""
            ID = ""
            SY = ""
            DI = ""
            atcclink = ""
            dsmzlink = ""
            doublingt = ""
            hasebv = False
            origin = ""
            isMeta = False
            characteristics = ""
            ancestry = ""
            issue = ""
            instability = ""
            transfected = ""
            depmap = set()
            lid = None
            individual = "ID-" + randomString(
                stringLength=6, stype="all", withdigits=True
            )
            comments = [""]
        elif line[:2] == "HI":
            parent = line[5:].split(" !")[0]
            if parent in cl:
                individual = cl[parent][5]
            else:
                tofind.update({lid: parent})
        elif line[:2] == "OX":
            if "Homo sapiens" in line:
                Homo = True
        elif line[:2] == "OI":
            name = line[5:].split(" !")[0]
            # adding its same id
            if name in cl:
                individual = cl[name][5]
            else:
                tofind.update({lid: name})
        elif line[:2] == "ID":
            ID = line[5:]
        elif line[:2] == "SY":
            SY = line[5:]
        elif line[:2] == "DI":
            DI += line[5:] + "; "
        elif line[:2] == "DT":
            date = line.split("Created: ")[-1].split(";")[0]
        elif line[:2] == "DR":
            if "DepMap; " in line:
                depmap.add(line.split("DepMap; ")[-1])
            if "ATCC" in line:
                atcclink = line.split("; ")[-1]
            if "DSMZ" in line:
                dsmzlink = line.split("; ")[-1]
        elif line[:2] == "CA":
            if "Cancer cell line" in line:
                CL = True
        elif line[:2] == "AG":
            age = line[5:]
            i = age.split("Y")[0]
            if i.endswith("FW"):
                # foetus
                i = "Fetus"
            elif i == "C":
                i = "Children"
            elif i in [
                "Embryonic stage",
                "Late embryonic stage",
                "Blastocyst stage",
                "Neonate larva",
            ]:
                i = "Embryo"
            elif i == "Age unspecified":
                i = "U"
            elif i not in ["Adult", "Children", "Fetus", "Embryo"]:
                if i.endswith("M") or i.endswith("W") or i.endswith("D"):
                    # 0Y
                    i = 0
                elif i.endswith("Y"):
                    i = i.split("-")[0].split("<")[-1].split(">")[-1]
                    i = float(i)
            age = i
        elif line[:2] == "SX":
            sex = line[5:]
            if sex == "":
                sex = "U"
            if sex == "Sex unspecified":
                sex = "U"
        elif line[:2] == "CC":
            comments.append(line[5:])
            if "Problematic cell line:" in line:
                issue += line.split("Problematic cell line: ")[-1] + " "
            if "Discontinued: DepMap; " in line:
                issue += "Removed from DepMap; "
                depmap.remove(line.split(";")[1][1:])
            if "Doubling time: " in line:
                doublingt = line.split("Doubling time: ")[-1].split("(")[0]
                if "everal week" in doublingt:
                    t = 13 * 24
                else:
                    t = float(
                        doublingt.split(" ")[0]
                        .split("-")[0]
                        .split("~")[-1]
                        .split(">")[-1]
                        .split("=")[-1]
                        .split("<")[-1]
                    )
                    if "day" in doublingt:
                        t = t * 24
                doublingt = int(t)
            if "Epstein-Barr virus (EBV)" in line:
                hasebv = True
            if "Derived from sampling site: " in line:
                origin = line.split("Derived from sampling site: ")[-1][:-1]
            if "Derived from metastatic site:" in line:
                origin = line.split("Derived from metastatic site: ")[-1][:-1]
                isMeta = True
            if "Characteristics: " in line:
                characteristics = line.split("Characteristics: ")[-1][:-1]
            if "Genome ancestry: " in line:
                ancestry = line.split("Genome ancestry: ")[-1][:-1]

            if "Caution: " in line:
                issue = line.split("Caution: ")[-1][:-1]
            if "Microsatellite instability: " in line:
                instability = line.split(" (")[1].split(")")[0]
            if "Transfected with: " in line:
                transfected += line.split("Transfected with: ")[-1][:-1] + " | "
    l.close()
    for k, v in tofind.items():
        if k in cl and v in cl:
            cl[k][5] = cl[v][5]
    cld = pd.DataFrame(
        data=cl.values(),
        columns=[
            "depmap_id",
            "id",
            "disease",
            "age",
            "sex",
            "patient_id",
            "parent_id",
            "date",
            "synonyms",
            "comments",
            "atcclink",
            "dsmzlink",
            "doublingt",
            "hasebv",
            "origin",
            "isMeta",
            "characteristics",
            "ancestry",
            "issue",
            "instability",
            "transfected",
        ],
        index=cl.keys(),
    )
    return cld


def getAncestry(
    df,
    col="ancestry",
    ancestry={
        "African": [],
        "Native American": [],
        "East Asian, North": [],
        "East Asian, South": [],
        "South Asian": [],
        "European, North": [],
        "European, South": [],
    },
):
    """getAncestry returns a df with the ancestries

    Args:
        df ([type]): [description]
        col (str, optional): [description]. Defaults to "ancestry".
        ancestry (dict, optional): [description]. Defaults to {'African':[], 'Native American':[], 'East Asian, North':[], 'East Asian, South': [], "South Asian":[], "European, North":[], 'European, South':[]}.

    Returns:
        pd.DataFrame: [description]
    """
    name = []
    for i, anc in df[col].iteritems():
        if anc == "":
            continue
        for place in anc.split("; "):
            ori, perc = place.split("%")[0].split("=")
            ancestry[ori].append(float(perc))
        name.append(i)
    return pd.DataFrame(data=ancestry, index=name)


def _fetchFromServer(ensemble_server, attributes):
    server = BiomartServer(ensemble_server)
    ensmbl = server.datasets["hsapiens_gene_ensembl"]
    res = pd.read_csv(
        io.StringIO(
            ensmbl.search({"attributes": attributes}, header=1).content.decode()
        ),
        sep="\t",
    )
    return res


def generateGeneNames(
    ensemble_server="http://nov2020.archive.ensembl.org/biomart",
    useCache=False,
    cache_folder="/".join(__file__.split("/")[:-3]) + "/",
    attributes=[],
):
    """generate a genelist dataframe from ensembl's biomart

    Args:
        ensemble_server ([type], optional): [description]. Defaults to ENSEMBL_SERVER_V.
        useCache (bool, optional): [description]. Defaults to False.
        cache_folder ([type], optional): [description]. Defaults to CACHE_PATH.

    Raises:
        ValueError: [description]

    Returns:
        [type]: [description]
    """
    attr = [
        "ensembl_gene_id",
        "clone_based_ensembl_gene",
        "hgnc_symbol",
        "gene_biotype",
        "entrezgene_id",
    ]
    assert cache_folder[-1] == "/"

    cache_folder = os.path.expanduser(cache_folder)
    createFoldersFor(cache_folder)
    cachefile = os.path.join(cache_folder, ".biomart.csv")
    if useCache & os.path.isfile(cachefile):
        print("fetching gene names from biomart cache")
        res = pd.read_csv(cachefile)
    else:
        print("downloading gene names from biomart")
        res = _fetchFromServer(ensemble_server, attr + attributes)
        res.to_csv(cachefile, index=False)

    res.columns = attr + attributes
    if type(res) is not type(pd.DataFrame()):
        raise ValueError("should be a dataframe")
    res = res[~(res["clone_based_ensembl_gene"].isna() & res["hgnc_symbol"].isna())]
    res.loc[res[res.hgnc_symbol.isna()].index, "hgnc_symbol"] = res[
        res.hgnc_symbol.isna()
    ]["clone_based_ensembl_gene"]

    return res


def cutLoops(li):
    """given a list of tuples representing edges of a graph like so: [(a,b),...], cut loops by removing tuples

    this list can represent potential graph like [(a,b),(b,c),(c,d),(c,a)], in this case this function would drop (c,a)
    """
    reses = []
    nli = []
    for j, (a, b) in enumerate(li):
        showcount(j, len(li))
        if a == b or (b, a) in nli:
            continue
        ini = False
        loc = -1
        drop = False
        for i, res in enumerate(reses):
            if a in res:
                if b in res:
                    # loop
                    drop = True
                    break
                else:
                    if ini:
                        # b was in another res
                        # merge b's group and a's group
                        poped = reses.pop(loc)
                        i = i - 1 if loc < i else i
                        reses[i] = reses[i] | poped
                        break
                    else:
                        ini = True
                        loc = i
                        reses[i].add(b)
            elif b in res:
                if ini:
                    # a was in another res, merge them
                    poped = reses.pop(loc)
                    i = i - 1 if loc < i else i
                    reses[i] = reses[i] | poped
                    break
                else:
                    ini = True
                    loc = i
                    reses[i].add(a)

        if not drop:
            nli.append((a, b))
            if not ini:
                reses.append(set([a, b]))
    return nli


def removeCoVar(mat, maxcorr=0.95):
    """removeCoVar list columns to remove as they covar with other columns

    just regular linear correlation.
    It displays a list of genes that have been dropped because their correlation
    to another gene is above a certain value.
    It shows a python dictionary {gene_to_be_dropped: gene_it_correlates_to}.

    Args:
        mat (array like): the matrix of obs x var
        maxcorr (float, optional): the max correlation above which to drop an observation. Defaults to 0.95.

    Returns:
        list(tuples): lists of observations to drop and their covarying observation to keep [(todrop,tokeep),...]
    """
    mat = mat.T
    loc = np.argwhere(np.corrcoef(mat) >= maxcorr)
    nloc = cutLoops(loc)

    drop = []
    sameness = []
    for a, b in nloc:
        if a not in drop:
            # we already said to drop b: do nothing here
            if b not in drop:
                drop.append(b)
                sameness.append((b, a))
    if type(mat) is pd.DataFrame:
        col = mat.index.tolist()
        # replace sameness values with the col values
        res = []
        for (i, j) in sameness:
            res.append((col[i], col[j]))
        sameness = res

    return sameness
