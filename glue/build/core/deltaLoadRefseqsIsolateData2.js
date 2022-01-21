// Load EVE data from tab file 
var loadResult;
glue.inMode("module/deltaTabularUtility", function() {
	loadResult = glue.tableToObjects(glue.command(["load-tabular", "tabular/core/deltaretrovirus-reference-data.tsv"]));
	//glue.log("INFO", "load result was:", loadResult);
});


_.each(loadResult, function(refSeqObj) {

    //glue.log("INFO", "refSeqObj was:", refSeqObj);
	glue.inMode("custom-table-row/isolate_data/"+refSeqObj.sequenceID, function() {
	
		glue.command(["set", "field", "host_sci_name", refSeqObj.host_sci_name]);
		glue.command(["set", "field", "host_common_name", refSeqObj.host_common_name]);

	});

	glue.inMode("sequence/ncbi-refseqs/"+refSeqObj.sequenceID, function() {
	
		glue.command(["set", "field", "name", refSeqObj.name]);
		glue.command(["set", "field", "full_name", refSeqObj.full_name]);
		glue.command(["set", "field", "subgenus", refSeqObj.subgenus]);
		glue.command(["set", "field", "clade", refSeqObj.clade]);
	});


});
