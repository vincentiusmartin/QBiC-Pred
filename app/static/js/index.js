function listPreddir(){
    // https://jsfiddle.net/qdmwxkb5/
    // https://stackoverflow.com/questions/815103/jquery-best-practice-to-populate-drop-down
    var $predDropdown = $("#predlist");
    var $famDropdown = $("#familylist");
    $.ajax({
      type: "GET",
      url: '/predlist',
      dataType: 'json',
      success: function(data){
          // parse the resulting json
          $.each(data, function(family, gene_pbm) {
            // Fill #predlist
              var $group = $('<optgroup label="' + family + '" />');
              $.each(gene_pbm, function(gene, pbmnames) {
                  // should be searchable by family and by TF-name
                  $group.append($('<option />').val(gene+":"+pbmnames).text(gene).attr("data-tokens",family+","+gene));  // family,gene
              })
              $predDropdown.append($group);

              // Fill #familylist
              $famDropdown.append($('<option />').val(family).text(family));
          });
          $('.selectpicker').selectpicker('refresh'); // needed to refresh optionlist
      },
      error: function() {
          alert("error getting prediction list");
      }
    });
}

function updateFromFamilies(){
    /* Token: {family:genename}*/
    /* When family list is closed */
    $('#familylist').on('hide.bs.select', function () {
        var selected = new Array();
        $("#familylist option:selected").each(function(){
            selected.push(this.value);
        });

        $('#predlist').selectpicker('deselectAll');
        toselect = [];
        for(var i = 0; i < selected.length; i++) {
            $('#predlist option').each(function() {
                if($(this).attr("data-tokens").split(",")[0] == selected[i]) {
                    toselect.push($(this).val());
                }
            });
        }
        $('#predlist').selectpicker('val',toselect);
    });

    /* When pred list is closed */
    $('#predlist').on('hide.bs.select', function () {
        var famSelectedArr = new Array();
        $("#predlist option:selected").each(function(){
            fam = $(this).attr("data-tokens").split(",")[0];
            if(!famSelectedArr.includes(fam)){
                famSelectedArr.push(fam);
            }
        });
        $('#familylist').selectpicker('deselectAll');
        $('#familylist').selectpicker('val',famSelectedArr);
    });
}

function updateOutlabel(){
    $("#outputradio").change(function() {
        var val = $('input[name=optradio]:checked', '#outputradio').val();
        var newOpts = [];
        if (val=="1"){
            $("#output-option").html("#TFs to output:");
            newOpts = [1,2,3];
        }else{
            $("#output-option").html("p-val cutoff:");
            newOpts = [0.01,0.05,0.1];
        }
        var $opt = $("#output-selection-opt");
        $opt.empty();
        $.each(newOpts,function(idx,val){
            $opt.append($("<option />").text(val));
        });
    });
}

/* get TF uploaded from file */
function uploadTFFomFile(){
    var file = $('#tf-file').prop('files')[0];
    if (file) {
        var reader = new FileReader();
        reader.readAsText(file, "UTF-8");
        reader.onload = function (evt) {
            var lines = evt.target.result.trim().split('\n');
            var selected = [];
            for(var i = 0; i < lines.length; i++) {
                var tfval = $('#predlist option').filter(function() {
                    return (this.text.localeCompare(lines[i]) == 0);
                }).val(); // ugly jquery syntax...
                selected.push(tfval);
            }
            $('#predlist').selectpicker('val',selected);
        }
        reader.onerror = function (evt) {
            alert("error");
        }
    }
    $('#tf-file').val('');
}

// --------- Functionality related ---------


function uploadFile(){
    //$('#upload-form').html("<p>checking file...</p><br />");
    var form = $('#input-form')[0];
    var fd = new FormData(form);

    $('#submit-job').prop('disabled', true).html("<i class=\"fas fa-circle-notch fa-spin\"></i> Submit");
    //$('#upload-msg').html("Uploading your file...").css("color","blue");
    $.ajax({
        type: 'POST',
        url: '/upload',
        data: fd, // pass the form data to Flask
        processData: false,  // tell jQuery not to process the data
        contentType: false,   // tell jQuery not to set contentType
        success: function(data, status, request) { // jsonified return from handle_upload in views.py
            // upload is successful
            $('#upload-msg').html(""); // empty error message
            status_url = request.getResponseHeader('Location');
            window.location.replace(status_url); // go to the result page, pass the url
        },
        error: function(data,status) {
            $('#submit-job').prop('disabled', false).html("Submit");
            //$('#upload-msg').show();
            $('#upload-msg').html(data["responseJSON"]["Message"]);
        }
    });
}

function testing(){
  $.ajax({
      type: "GET",
      url: '/testing',
      dataType: 'json',
      success: function(res){
            alert(res);
            var tbl = '<table class=\"table table-bordered\" id=\"testtbl\">';
            var keys = [];
            tbl += '<thead>';
            tbl += '<tr>';
            for (var i = 0; i < Object.keys(res[0]).length; i++) { // zero index is the list of colname
                keys.push(res[0][i]);
                tbl += '<th>' + res[0][i] + '</th>';
            }
            tbl += '</tr></thead><tbody>';
            for (var i = 1; i < res.length; i++) {
                tbl += '<tr>';
                for (var j = 0; j < keys.length; j++){
                    tbl += '<td>' + res[i][keys[j]] + '</td>';
                }
                tbl += '</tr>';
            }
            tbl += '</tbody></table>';
            $('#testing').html(tbl);
            //$('#testtbl').DataTable();
      },
      error: function() {
          alert("error testing");
      }
    });
}

// --------- End of Progress Related ---------

// jquery specific function
$(function() {
    listPreddir(); // list all TFs from our list in the dropdowns
    updateOutlabel(); // check change on output options
    updateFromFamilies();
    $('#submit-job').click(uploadFile);
    $('#submit-tf').click(uploadTFFomFile);
    //testing();

    $('[data-toggle="popover"]').popover();
});
