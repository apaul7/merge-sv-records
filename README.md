#sv-merge
The purpose of this tool is to locate deletion and duplication calls that have been split up by the caller and recover them. The idea is that a window size can be used so that if sv1 end is within a window of sv2 start the event can be merged a new record put in the output vcf. 

Future?
- add tests
- confirm no other calls within window
- redo depth coverage
- 
