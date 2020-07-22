// Load EVE data from tab file 
var loadResult;
glue.inMode("module/deltaTabularUtility", function() {
	loadResult = glue.tableToObjects(glue.command(["load-tabular", "tabular/erv/erv-locus-data.tsv"]));
	//glue.log("INFO", "load result was:", loadResult);
});

_.each(loadResult, function(ervObj) {

	glue.inMode("custom-table-row/locus_data/"+ervObj.sequenceID, function() {
	
		glue.command(["set", "field", "scaffold", ervObj.scaffold]);
		glue.command(["set", "field", "start", ervObj.extract_start]);
		glue.command(["set", "field", "end", ervObj.extract_end]);
		glue.command(["set", "field", "orientation", ervObj.orientation]);

	});

	glue.inMode("sequence/fasta-digs-erv/"+ervObj.sequenceID, function() {
	
		glue.command(["set", "field", "name", ervObj.sequenceID]);
		glue.command(["set", "field", "full_name", ervObj.sequenceID]);
		glue.command(["set", "field", "species", ervObj.organism]);

		glue.command(["set", "field", "length", ervObj.sequence_length]);
	});
	

});
