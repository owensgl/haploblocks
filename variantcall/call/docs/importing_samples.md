[Back to start](./README.md)


Importing Samples
====================

Most of the samples we've used have been sequenced on NanuQ, which
produces FastQ files and their associated MD5s. Consequently, the
instructions here will be Nanuq centric.

Here is an example sample definition specifying the location of FastQ files and their MD5s.

```yaml
lane: HI.4530

samples:
  - name: ANN0803
    locations:
      fastq:
        - https://genomequebec.mcgill.ca/nanuqMPS/fileDownload/id/417539/type/READ_SET_FASTQ/filename/HI.4530.004.index_22.ANN0803_R1.fastq.gz
        - https://genomequebec.mcgill.ca/nanuqMPS/fileDownload/id/417539/type/READ_SET_FASTQ_PE/filename/HI.4530.004.index_22.ANN0803_R2.fastq.gz
      md5:
        - https://genomequebec.mcgill.ca/nanuqMPS/readSetMd5Download/id/417539/type/READ_SET_FASTQ/filename/HI.4530.004.index_22.ANN0803_R1.fastq.gz.md5
        - https://genomequebec.mcgill.ca/nanuqMPS/readSetMd5Download/id/417539/type/READ_SET_FASTQ_PE/filename/HI.4530.004.index_22.ANN0803_R2.fastq.gz.md5
```

You will notice that this file is somewhat tedious to specify, so in
the repository there are tools which will convert a download list from
NanuQ into usable sample definitions.

Getting the list of URLs
---------------

To import samples from Nanuq, what you need to know is:
 1. the username and password to the genomequebec [web interface](https://genomequebec.mcgill.ca/nanuqAdministration/nanuq-administration/welcomeHome.do)
 2. the readset links zip file.

To get the latter (2):

 1. Login to Nanuq. Pick your project.
 1. Click the "HiSeq Read Sets" tab.
 1. Click the top checkbox to mark ALL elements. Increase the number shown per page if necessary so you can see them all at once.
 1. Click "Download Read Files" button at the top.
 1. Select "Download files from selected reads" with the radiobutton in the popup.
 1. Choose "Compressed file", and select all the checkboxes at the bottom.
 1. Click "Download". save the zip file somewhere.

Inside that zip file is a file called "readSetLinks.txt". That's the file you want.

It will contain one URL per line:
```
https://genomequebec.mcgill.ca/nanuqMPS/fileDownload/id/450118/type/READ_SET_FASTQ_PE/filename/HI.4663.002.Rieseberg_2.DEB_896_R2.fastq.gz
https://genomequebec.mcgill.ca/nanuqMPS/fileDownload/id/450118/type/READ_SET_FASTQ/filename/HI.4663.002.Rieseberg_2.DEB_896_R1.fastq.gz
https://genomequebec.mcgill.ca/nanuqMPS/fileDownload/id/450119/type/READ_SET_FASTQ_PE/filename/HI.4663.002.Index_27.DEB_1837_R2.fastq.gz
https://genomequebec.mcgill.ca/nanuqMPS/fileDownload/id/450119/type/READ_SET_FASTQ/filename/HI.4663.002.Index_27.DEB_1837_R1.fastq.gz
...
https://genomequebec.mcgill.ca/nanuqMPS/readSetMd5Download/id/450120/type/READ_SET_FASTQ_PE/filename/HI.4663.002.Index_15.664647_GIG_R2.fastq.gz.md5
...
```

Converting the URL list into usable sample definitions.
---------------------------------------------------

It is tedious to work with these URL files. To convert that list of
FastQ and MD5 URLs into a usable sample definition file, we run the
command:

```bash

# make sure you have the genomics conda environment activated.

cd <repo>/snake/
./scripts/process-sample-list.py readSetLinks.txt --mode yaml > new-samples.yaml

```

You may inspect the produced sample file, or annotate it with additional information.

 1. You should make sure you are not creating duplicates. This may happen if the new list of samples is an addition to a previously generated list. (e.g. as new samples arrive in batches).
 2. If your sample names have special URL characters in them, such as `?`, or `#`, you should
    modify the sample name information (but leave the URL as it was given to you).
    
When you are ready to use it, place the new file in the `<repo>/snake/samples/`
folder.

Pre-downloading the files in the listing
-------------------------------------------------

Most likely, you will be running multiple alignments in parallel, so
you will benefit from pre-downloading sample files before executing
the pipeline. Nanuq, in particular, can only support up to 10
simultaneous downloads, realistically. If you ask for more, you will
get errors.

We can use two other scripts for this.

 1. convert the readSetLinks into a format that the downloader program will accept.

    ```bash
    ./scripts/process-sample-list.py readSetLinks.txt --mode urls > new-urls.txt
    ```

 1. Run the downloader script to read that `new-urls.txt` file, and
 download to a cache folder. You want to specify the same cache folder
 that is defined inside your `config.yaml`, under key `casdir`.

    ```bash
    
    cd <repo>/snake

    # check the CASDIR you are using in the pipeline
    cat config.yaml | grep casdir
    
    # replace CASDIR with the directory found above
    ./scripts/download-samples.py < ./scripts/new-urls.txt --creds creds.yaml data/nanuq_dl/ --workdir data/tmp/
    ```
    
    In the cache folder, files are indexed by their MD5, so that the
    bulk entire downloads can be skipped if the content is already
    known to the directory. So, don't worry if you've partly
    downloaded the files already, those downloads will be skipped.

    All tools from the pipelines which necessitate downloading samples
    can optionally take in a cache directory created earlier.
