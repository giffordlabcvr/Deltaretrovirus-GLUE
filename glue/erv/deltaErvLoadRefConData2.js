// Load ERV data from tab file 
var loadResult;
glue.inMode("module/deltaTabularUtility", function() {
	loadResult = glue.tableToObjects(glue.command(["load-tabular", "tabular/erv/erv-refseqs-side-data.tsv"]));
	// glue.log("INFO", "load result was:", loadResult);
});

_.each(loadResult, function(ervObj) {

	glue.inMode("custom-table-row/refcon_data/"+ervObj.sequenceID, function() {
	
		glue.log("INFO", "Entering locus data for ERV reference:", ervObj.sequenceID);
		glue.command(["set", "field", "reftype", ervObj.reftype]);

	});

	glue.inMode("sequence/fasta-refseqs-erv/"+ervObj.sequenceID, function() {
	
		glue.log("INFO", "Entering sequence table data for ERV reference:", ervObj.sequenceID);

		glue.command(["set", "field", "name", ervObj.name]);
		glue.command(["set", "field", "full_name", ervObj.full_name]);
		glue.command(["set", "field", "subgenus", ervObj.subgenus]);
		glue.command(["set", "field", "host_group_name", ervObj.host_group_name]);
		glue.command(["set", "field", "host_group_taxlevel", ervObj.host_group_taxlevel]);


	});

});


