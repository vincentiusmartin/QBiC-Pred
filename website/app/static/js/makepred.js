function changeInputOptOnClick(){
    $(".disabled-btn").click(function(){
        if ($(".disabled-btn").val() == 1) {
            changeInputOpt(1);
        }else{
            changeInputOpt(2);
        }
    });
}

function changeInputDownloadLink(){
    var filename = $("#examplelist option:selected").attr("inputfile");
    if(filename){
        $("#download-input-select").attr("href","/download/"+filename).addClass("grey-download").html("Download mutation dataset");
    }else{
        $("#download-input-select").attr("href","").html("");
    }
}

function changeOptOnFileSelected(){
    /* Also handle for file selection here */
    $("#input-file").on('change', function(){
        if($("#input-mode").val() == "2"){
            changeInputOpt(1);
        }
    });
}

/* Make separate function with on click so we can independently call it */
function changeInputOpt(opt){
    var content1 = "Upload a mutation dataset:";
    var content2 = "Select an example use case <br /> (this will fill all necessary fields): ";
    if(opt == 2){ // 2 is active and 1 is inactive
        $("#input-mode").val("2"); // needed to differentiate input type
        $("#input-lbl-2").html(content2);
        $("#input-lbl-1").html('<button class="btn disabled-btn" value="1" type="button">' +
              content1 +
            '</button>');
        var optSelected = $("#examplelist option:selected");
        changeInputDownloadLink();
        $("#examplelist").change(function(){
            changeInputDownloadLink();
        });
    }else{
        $("#input-mode").val("1");
        $("#input-lbl-1").html(content1);
        $("#input-lbl-2").html('<button class="btn disabled-btn" value="2" type="button">' +
              content2 +
            '</button>');
        $("#download-input-select").attr("href","").html("");
    }
    changeInputOptOnClick(); // reload java script
}

function listExampleInput(){
    var dfd = jQuery.Deferred();

    var $exampleDropdown = $("#examplelist");
    /* Fill the example dropdown */
    var $noneOpt = $('<option />').val("").text("None").attr("inputfile","")
              .attr("tfs",[]).attr("genomever","hg19").attr("outputtype",1);
    $exampleDropdown.append($noneOpt); // none option
    $.ajax({
      type: "GET",
      url: '/examplelist',
      dataType: 'json',
      success: function(data){
          dfd.resolve("ok");
          // parse the resulting json
          $.each(data, function(name,attributes) {
              var $group = $('<option />').val(name).text(name);
              $.each(attributes, function(key,val) {
                  $group.attr(key,val);
              });
              $exampleDropdown.append($group);
          });
          $exampleDropdown.selectpicker('refresh');
      },
      error: function() {
          dfd.reject("error");
          alert("error getting example list");
      }
    });

    /* add action when an option is selected */
    // predlist
    $exampleDropdown.on("changed.bs.select",function(){
        var optSelected = $("#examplelist option:selected");
        /* Update tf list */
        var tfs = optSelected.attr("tfs").split(",");
        var tfsVal = []
        for(var i = 0; i < tfs.length; i++) {
            // get value from its text
            var curVal = $('#predlist option').filter(function () { return $(this).html() == tfs[i]; }).val();
            tfsVal.push(curVal);
        }
        $('#predlist').selectpicker('deselectAll');
        $('#predlist').selectpicker('val',tfsVal);
        updateToFamilies();

        /* Update genome version */
        $('#genomelist').selectpicker('val',optSelected.attr("genomever"));
        /* Update output type */
        $('input[name=optradio][value=' + optSelected.attr("outputtype") + ']').prop('checked',true)
        updateOutputLabel(optSelected.attr("outputtype"));

        /* Update input */
        if(optSelected.val().length > 0){
            changeInputOpt(2);
        }else {
            changeInputOpt(1);
        }
    });
    return dfd.promise();
}

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
          $famDropdown.selectpicker('refresh'); // needed to refresh optionlist
          $predDropdown.selectpicker('refresh');
      },
      error: function() {
          alert("error getting prediction list");
      }
    });
}

function makeTFDownloadLink(){
    $("#predlist").change(function(){
        var tfSelectedArr = new Array();
        $("#predlist option:selected").each(function(){
            tfSelectedArr.push($(this).text());
        });
        if(tfSelectedArr.length>0){
            $("#download-tfs-select").attr("href","/tfsdownload/"+JSON.stringify(tfSelectedArr)).addClass("grey-download").html("Download selected TFs");
        }else{
            $("#download-tfs-select").attr("href","").html("");
        }
    });
}

function updateToFamilies(){
    var famSelectedArr = new Array();
    $("#predlist option:selected").each(function(){
        fam = $(this).attr("data-tokens").split(",")[0];
        if(!famSelectedArr.includes(fam)){
            famSelectedArr.push(fam);
        }
    });
    $('#familylist').selectpicker('deselectAll');
    $('#familylist').selectpicker('val',famSelectedArr);
}

function updateFromFamilies(){
    /* Token: {family:genename}*/
    /* When family list is closed */
    $('#familylist').on('hide.bs.select', function () { //hide.bs.select
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
        updateToFamilies()
    });
}

function updateOutputLabelWrapper(){
    $("#outputradio").change(function() {
        var val = $('input[name=optradio]:checked', '#outputradio').val();
        updateOutputLabel(val);
    });
}

function updateOutputLabel(newlabel){
    if (newlabel=="1"){
        $("#output-option").html("#TFs to output:");
        $("#output-selection-wrapper").html(`
            <select id="output-selection-opt" name="output-selection-opt">
              <option>1</option>
              <option>2</option>
              <option>3</option>
            </select>
            `)
    }else{
        $("#output-option").html("p-value cutoff:");
        $("#output-selection-wrapper").html(`
            <input class="form-control" type="text"
             id="output-selection-opt" name="output-selection-opt"
             value="0.0001" placeholder="input p-value cutoff" />
            `).addClass("col-lg-6");
    }
}

/* get TF uploaded from file */
function uploadTFFomFile(){
    $("#tf-file").on('change', function(){
        var file = $('#tf-file').prop('files')[0];
        if (file) {
            var reader = new FileReader();
            reader.readAsText(file, "UTF-8");
            reader.onload = function (evt) {
                var lines = evt.target.result.trim().split('\n');
                $.ajax({
                    type: 'POST',
                    url: '/checktfnames',
                    contentType: 'application/json',
                    dataType : 'json',
                    data: JSON.stringify({'tfs':lines}),
                    success: function(data,status,request){
                        // parse the resulting json
                        var nothgnc = request.getResponseHeader('notfound').split(",");
                        var inhgnc = request.getResponseHeader('found').split(",");
                        var notAvail = [];
                        var selected = [];
                        for(var i = 0; i < inhgnc.length; i++) {
                            var tfname = inhgnc[i].trim();
                            var tfval = $('#predlist option').filter(function() {
                                return (this.text.localeCompare(tfname) == 0);
                            }).val(); // ugly jquery syntax...
                            if (typeof tfval === "undefined") {
                                notAvail.push(inhgnc[i]);
                            }else{
                                selected.push(tfval);
                            }
                        }
                        if(notAvail.length > 0 || nothgnc.length > 0){
                            alert("Not all TFs could be processed:\n" +
                            selected.length + " input TFs have PBM data and were selected.\n" +
                            notAvail.length + " input TFs do not have PBM data:" + notAvail.toString() + "\n" +
                            nothgnc.length + " input TFs are probably not using the HGNC nomenclature for gene names: " +
                            nothgnc.toString()+
                            "\n(TF gene names must use HGNC nomenclature)");
                        }
                        $('#predlist').selectpicker('val',selected);
                        updateToFamilies();
                    },
                    error: function() {
                        alert("error in reading the TF list");
                    }
                });
            }
            reader.onerror = function (evt) {
                alert("file error");
            }
        }
    });
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

// --------- End of Progress Related ---------

// jquery specific function
$(function() {
    $("#input-file").filestyle({
        btnClass: 'btn-secondary',
        htmlIcon: '<i class="fas fa-folder-open"></i> ',
        size:'sm'
    });

    listPreddir(); // list all TFs from our list in the dropdowns
    makeTFDownloadLink();
    changeOptOnFileSelected();
    listExampleInput();
    changeInputOptOnClick();

    updateOutputLabelWrapper(); // check change on output options
    updateFromFamilies();
    uploadTFFomFile();
    $('[data-toggle="popover"]').popover();

    // submit the form
    $('#submit-job').click(uploadFile);

    $('.navbar-nav .nav-link').removeClass('active');
    $('#nav-2').addClass('active');
});
