# Build Genome Datalake Tables

I want to add a new function to the spec and UI called "build_genome_datalake_tables". This function should take as input "input_refs", which is a list of workspace references references to Genome or GenomeSet. The output should be a object:
typedef structure {
        string report_name;
        string report_ref;
        string workspace;
    } ReportResults;
The object should be named "output". The UI should accept a list of Genome or GenomeSet workspace objects. A multi-input should be used. We should also have as input: a string suffix, a boolean indicating if the generated models should be saved.
