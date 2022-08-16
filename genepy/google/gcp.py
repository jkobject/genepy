# GCPFunction.py

from google.cloud import storage
import os
import subprocess
import re
from genepy.utils import helper as h
import signal


def list_blobs_with_prefix(bucket_name, prefix, delimiter=None):
    """Lists all the blobs in the bucket that begin with the prefix.

    This can be used to list all blobs in a "folder", e.g. "public/".

    The delimiter argument can be used to restrict the results to only the
    "files" in the given "folder". Without the delimiter, the entire tree under
    the prefix is returned. For example, given these blobs:

            /a/1.txt
            /a/b/2.txt

    If you just specify prefix = '/a', you'll get back:

            /a/1.txt
            /a/b/2.txt

    However, if you specify prefix='/a' and delimiter='/', you'll get back:

            /a/1.txt

    """
    storage_client = storage.Client()
    bucket = storage_client.get_bucket(bucket_name)
    ret = []
    blobs = bucket.list_blobs(prefix=prefix, delimiter=delimiter)
    for blob in blobs:
        ret.append(blob.name)
    return ret


def mvFiles(files, location, group=50, listen_to_errors=False):
    """
    move a set of files in parallel (when the set is huge)

    Args:
    ----
        files: gs paths
        location: to move the files to
        group: files to do in parallel
    """
    by = len(files) if len(files) < group else group
    for sfiles in h.grouped(files, by):
        a = ""
        for val in sfiles:
            a += val + " "
        code = os.system("gsutil -m mv " + a + location)
        if code != 0 and listen_to_errors:
            print("pressed ctrl+c or command failed")
            break


def lsFiles(files, add="", group=50):
    """
    list a set of files in parallel (when the set is huge)

    Args:
    ----
        files: gs paths
        add: additional params to add
        group: files to do in parallel
    """
    print("listing files in gs")
    by = len(files) if len(files) < group else group
    res = []
    for sfiles in h.grouped(files, by):
        a = ""
        for val in sfiles:
            a += val + " "
        data = subprocess.run(
            "gsutil -m ls " + add + " '" + a + "'", capture_output=True, shell=True
        )
        if data.returncode != 0:
            if "One or more URLs matched no objects" not in str(data.stderr):
                raise ValueError("issue with the command: " + str(data.stderr))
        if len(str(data.stdout)) < 4:
            return []
        res += (
            str(data.stdout)[2:-1].split("\\n")[:-1]
            if "L" not in add
            else ["gs://" + i for i in str(data.stdout).split("\\ngs://")]
        )
        if "TOTAL:" in res[-1] and "L" not in add:
            res = res[:-1]
    return res


def cpFiles(files, location, group=50):
    """
    copy a set of files in parallel (when the set is huge)

    Args:
    ----
        files: gs paths
        location to copy
        group: files to do in parallel
    """
    by = len(files) if len(files) < group else group
    for sfiles in h.grouped(files, by):
        a = ""
        for val in sfiles:
            a += val + " "
        code = os.system("gsutil -m cp " + a + location)
        if code != 0:
            print("pressed ctrl+c or command failed")
            break


def catFiles(files, group=50, split=False, cut=False):
    """
    copy a set of files in parallel (when the set is huge)

    Args:
    ----
        files: gs paths
        location to copy
        group: files to do in parallel
        cut: split all lines into chunks of size cut
        split: split lines by split e.g. \\n
    """
    by = len(files) if len(files) < group else group
    res = []
    for i, sfiles in enumerate(h.grouped(files, by)):
        print(i / (len(files) / by))
        a = ""
        for val in sfiles:
            a += val + " "
        data = subprocess.run("gsutil -m cat " + a, capture_output=True, shell=True)
        if data.returncode != 0:
            if "One or more URLs matched no objects" not in str(data.stderr):
                print(ValueError("issue with the command: " + str(data.stderr)))
                return res
        if len(str(data.stdout)) < 4:
            return []
        resa = str(data.stdout)[2:-1]
        if cut:
            res += [resa[i * cut : (i + 1) * cut] for i in range(int(len(resa) / cut))]
        elif split:
            res += resa.split(split)
        else:
            res += [resa]
    return res


def rmFiles(files, group=50, add="", dryrun=True):
    """
    remove a set of files in parallel (when the set is huge)

    Args:
    ----
        files: gs paths
        group: number to do in parallel
        add: additional gsutil cp params
    """
    by = len(files) if len(files) < group else group
    for sfiles in h.grouped(files, by):
        a = ""
        for val in sfiles:
            a += " " + val
        if add:
            add = " " + add
        if dryrun:
            print("gsutil -m rm" + add + a)
        else:
            code = os.system("gsutil -m rm" + add + a)
            if code != 0:
                print("pressed ctrl+c or command failed")
                break


def recoverFiles(files, cores=1):
    """
    recover a set of files in parallel that were erased

    files need to have their #id appended found using ls -al file

    Args:
    ----
        files: gs paths
        location: to move the files to
    """
    cmd = ["gsutil mv " + f + " " + f.split("#")[0] for f in files]
    h.parrun(cmd, cores=cores)


def folderRN(gspath, newpath, cores=1):
    """ """
    lis = lsFiles([gspath])
    if lis != 0:
        h.parrun(["gsutil -m mv " + val + " " + newpath for val in lis], cores=cores)
    else:
        raise ValueError("no such folder")


def patternRN(
    rename_dict,
    location,
    wildcards,
    types=[],
    dryrun=True,
    check_dependencies=True,
    cores=1,
):
    """
    rename/move a bunch of GCP objects found in some specific places

    Args:
    -----
        rename_dict: dict(prevName,newName)
        location:
        wildcards: list[str] can be one of  ['**', '.*', '*.','-.*'] if needs to be
                    ** means any occurence of this file in any folder will change its name
                    .* means all file unregarding of the suffix, will rename them all a.bam [a]da.bai to b.bam, [b]da.bai
                    *. means all files with the suffix, will change the suffix of these files from a to b
                    -.* means all file unregarding of the suffix, will rename them. not just replacing the a part with a to b but the full file name [a]dea.bam to b.bam
        types: Nothing yet
        test: if test, just shows the command but does not run it
        cores:  cores tells on how many processor to parallelize the tas#k
    """
    val = []
    for k, v in rename_dict.items():
        val.append(v)
        if k in val and check_dependencies:
            raise ValueError("circular dependency in the rename with key " + k)
    for k, v in rename_dict.items():
        loc = location
        if "**" in wildcards:
            loc += "**/"
        if "*." in wildcards or "-.*" in wildcards:
            loc += "*"
        loc += k
        if ".*" in wildcards or "-.*" in wildcards:
            loc += "*"
        res = os.popen("gsutil -m ls " + loc).read().split("\n")[:-1]
        print("found " + str(len(res)) + " files to rename")
        if "-.*" in wildcards:
            cmd = [
                "gsutil mv "
                + val
                + " "
                + "/".join(val.split("/")[:-1])
                + "/"
                + v
                + "."
                + ".".join(val.split("/")[-1].split(".")[1:])
                for val in res
            ]
        else:
            cmd = ["gsutil mv " + val + " " + val.replace(k, v) for val in res]
        if dryrun:
            print(cmd)
        else:
            h.parrun(cmd, cores=cores)


def get_all_sizes(folder, suffix="*"):
    """
    will sort and list all the files by their sizes.

    If some files have the same size, will list them together

    Args:
    ----
            folder: gs folder path
            suffix: of a specific file type

    Returns:
    -------
            dict(sizes:[paths])
    """
    samples = os.popen("gsutil -m ls -al " + folder + "**." + suffix).read().split("\n")
    # compute size filepath
    sizes = {
        "gs://"
        + val.split("gs://")[1].split("#")[0]: int(
            re.split("\d{4}-\d{2}-\d{2}", val)[0]
        )
        for val in samples[:-2]
    }
    names = {}
    for k, val in sizes.items():
        if val in names:
            names[val].append(k)
        else:
            names[val] = [k]
    if names == {}:
        # we didn't find any valid file paths
        print("We didn't find any valid file paths in folder: " + str(folder))
    return names


def exists(val):
    """
    tells if a gcp path exists
    """
    if type(val) is str:
        return os.popen("gsutil ls " + val).read().split("\n")[0] == val
    elif type(val) is list:
        rest = set(val) - set(lsFiles(val))
        return len(rest) == 0, rest


def extractSize(val):
    """
    extract the size from the string returned by an ls -l|a command
    """
    return "gs://" + val.split("gs://")[1].split("#")[0], int(
        re.split("\d{4}-\d{2}-\d{2}", val)[0]
    )


def extractTime(val):
    """
    extract the size from the string returned by an ls -l|a command
    """
    return val.split("  ")[1].split("T")[0]


def extractPath(val):
    """
    extract the path from the string returned by an ls -l|a command
    """
    return "gs://" + val.split("gs://")[1].split("#")[0]


def extractHash(val, typ="crc32c"):
    """
    extract the crc32 from the string returned by an ls -L command

    Args:
    ----
        type: flag ['crc32c','md5']
    """
    if "    Hash (crc32c):" in val and typ == "crc32c":
        return (
            val.split("    Hash (crc32c):          ")[-1]
            .split("\\\\n")[0]
            .split("\\n")[0]
        )
    elif "    Hash (md5):" in val and typ == "md5":
        return (
            val.split("    Hash (md5):          ")[-1].split("\\\\n")[0].split("\\n")[0]
        )
    else:
        return None


async def shareFiles(flist, users):
    """
    will share a list of files from gcp with a set of users.

    Args:
    ----
      users: list[str] of users' google accounts
      flist: list[str] of google storage path for which you want to share data

    """
    if type(users) is str:
        users = [users]
    for user in users:
        files = ""
        for i in flist:
            files += " " + i
        code = os.system("gsutil -m acl ch -ru " + user + ":R " + files)
        if code == signal.SIGINT:
            print("Awakened")
            break
    print("the files are stored here:\n\n")
    print(flist)
    print("\n\njust install and use gsutil to copy them")
    print("https://cloud.google.com/storage/docs/gsutil_install")
    print("https://cloud.google.com/storage/docs/gsutil/commands/cp")


def deleteOldVersions(path, onlymetagene=None, **kwargs):
    """ """
    data = subprocess.run("gsutil -m ls -alh " + path, capture_output=True, shell=True)
    if data.returncode != 0:
        if "One or more URLs matched no objects" not in str(data.stderr):
            print(ValueError("issue with the command: " + str(data.stderr)))
            return res
    if len(str(data.stdout)) < 4:
        return []
    resa = str(data.stdout)[2:-1].split("\\n")[:-2]
    torm = []
    for i, val in enumerate(resa):
        if onlymetagene is not None:
            if "metageneration=" + str(metagene) in val:
                name = "gs://" + val.split(" gs://")[1].split("  ")[0]
                torm.append(name)
            else:
                prevname = torm[-1].split("#")[0]
                name = "gs://" + val.split(" gs://")[1].split("#")[0]
                if prevname != name:
                    print(prevname + " is unique and won't be deleted")
                    torm.pop()
        else:
            if "metageneration=" + str(1) not in val:
                name = "gs://" + val.split(" gs://")[1].split("  ")[0]
                torm.append(name)
                if (
                    "gs://" + resa[i + 1].split(" gs://")[1].split("#")[0]
                    != name.split("#")[0]
                ):
                    print(name + " is unique and won't be deleted")
                    torm.pop()
    print(h.dups([val.split("#")[0] for val in torm]))

    return rmFiles(torm, **kwargs)
