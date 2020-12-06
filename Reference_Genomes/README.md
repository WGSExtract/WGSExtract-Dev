# Reference_Genomes subsystem

In releases v1 and v2, this was a static directory of delivered files.  And, in fact, was more than 90% of the release size content -- even though in BGZF compressed format.
Scripts in this directory may be incorporated into the Python code or remain as stand-alone utilities that the Python code and others can run. The Python code is being
modularized and versions for it to be stand-alone exist.

As it ends up, this whole system needed to be reworked. Especially the method of identifying reference genomes and then using them.  
This became critical as the needs to process CRAM files, not just BAMs, arose.  Especially after Nebula Genomics started delivering CRAM files.  But also because the
1K Genome project release samples and other ancient WGS sample sources were being delivered in CRAM format.

* [BAM Reference Model, Determing Your](https://bit.ly/34CO0vj) -- Google Docs document detailing the study (as it developed and was edited in)
* [Reference Genome Characteristics Spreadsheet](https://bit.ly/2ZmYPAg) -- Google Drive hosted LibreOffice Calc spreadsheet of reference genome characteristics 
(to be read in and used within the **WGS Extract** Python code)
* See the file **Process_Reference_Genomes.sh** file here for how we handle native downloaded files (as given in the spreadsheet link) to transform them into what is needed by the tools

This is still being worked into the actual Python code of the v3 release.  The results and files here are more a result of the detailed study into the issue and how to resolve
accurately identifying what a reference genome file is and contains; other than by an alluded file name reference that may or may not be captured in the BAM header.

Keys to note:
* No reference genome will be delivered with the v3 tools.  Cutting the release size from nearly 5 gigabytes down to 500 megabytes on this changealone
* The previoous reference genomes in the v2B release were inaccurately named, identified (typed) and used in the v2 code. Only one is so unique as it cannot be replicated
and must be retained in case someone used it outside the tool.
