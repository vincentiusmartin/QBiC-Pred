function displayOpt(searchOpt,searchKey){
    switch(searchOpt) {
        case 'exclude':
            return "~"+searchKey;
        case 'in sequence':
            return '&#x2731;'+searchKey+'&#x2731;';
        case 'at least':
            return '&geq;'+searchKey;
        case 'at most':
            return '&leq;'+searchKey;
        default:
            return searchKey;
    }
}

function addFilterElm(searchOpt,searchKey,searchCol){
    var ftext = searchCol + ":" + displayOpt(searchOpt,searchKey);
    $("#filter-list li").each(function(idx,li){
        if($(li).attr("searchKey")==searchKey && $(li).attr("searchCol")==searchCol){
            $(li).remove();
        }
    });
    var closeBtnText = "<span class=\"close\">x</span>"
    var $filterItem = $('<li />')
              .html(ftext + closeBtnText)
              .addClass("list-group-item")
              .addClass("list-group-item-custom")
              .attr("searchOpt",searchOpt)
              .attr("searchKey",searchKey)
              .attr("searchCol",searchCol);
    $("#filter-list").append($filterItem);
    attachCloseListHandler();
}

function getFilterList(){
    var filterList = [];
    $("#filter-list li").each(function(idx,li){
        var filter = {
            "searchOpt":$(li).attr("searchOpt"),
            "searchKey":$(li).attr("searchKey"),
            "searchCol":$(li).attr("searchCol")
        };
        filterList.push(filter);
    });
    return filterList;
}

function attachCloseListHandler(){
    $(".list-group-item").on( "click", ".close", function(){
        var parent = $(this).parent();
        var searchOpt = parent.attr("searchOpt");
        parent.remove();
        $('#restbl').DataTable().ajax.reload();
    });
}

/**
  This function parse result from datatable and display it to the user.
**/
function displayResult(status_url){
    var keys = [];
    // show result
    var job_arr = status_url.split("/");
    var job_id = job_arr[job_arr.length-1];
    var cols;
    $.getJSON('/getrescol/'+job_id,function(cols) {
        var searchFilter = [];
        $('#tablewrap').css('display','block');
        var restbl = $('#restbl').DataTable({
            dom: 'lr<"#ressearch">tip',
            // https://stackoverflow.com/questions/17237812/datatable-jquery-table-header-width-not-aligned-with-body-width
            //scrollX: true, -- don't use scrollX
            searching: false,
            processing: true, // Enable the display of 'processing' when the table is being processed
            serverSide: true,
            orderMulti: false,
            autoWidth: false,
            language: {
                processing: "<i class=\"fas fa-circle-notch fa-spin fa-5x\" style=\"color: gray;\"></i>"
            }, // $("#tableid").addClass("disabled");
            preDrawCallback: function( settings ) {
                /* disable some elments before draw */
                $('.dataTables_paginate').hide();
                $('.dataTables_length').hide();
                $('#ressearch-btn').prop('disabled', true);
            },
            drawCallback: function(settings, json) {
                /* reenable elments after draw */
                $('.dataTables_paginate').show();
                $('.dataTables_length').show();
                $('#ressearch-btn').prop('disabled', false);
                $(".cell-filter-item").click(function(){
                    searchOpt = $(this).data("filter");
                    searchKey = $(this).parent().siblings("button.cell-btn").text().trim();
                    searchCol = $(this).data("colname");
                    addFilterElm(searchOpt,searchKey,searchCol);
                    restbl.ajax.reload();
                });
            },
            ajax: {
                url: '/getrestbl/'+job_id,
                data:function(d){
                    searchFilter = getFilterList(); // get all filters from the list
                    d.searchFilter = JSON.stringify(searchFilter); // custom search
                }
            },
            columns: cols,
        });
        // this have to be outside due to the need of refering to the table
        $('#restbl').css({"width":"100%"}); // need to set this to align the header
        $('.dataTable').wrap('<div class="dataTables_scroll" />'); // this is needed to replace dataTables_scroll
        // This is used to replace the search bar from DataTable so we can have
        // our custom. Search bar itself is inactivated.
        // https://stackoverflow.com/questions/43454575/call-datatable-ajax-call-on-custom-button
        $("#ressearch").addClass("form-inline mb-3").html(`
          <div class="form-group" id="math-compare-form" style="margin-right:5px;display:none;">
            <select id="math-compare-select" name="math-compare-select" class="selectpicker"
              data-width="fit">
              <option selected>at most</option>
              <option>at least</option>
            </select>
          </div>
          <div class="form-group" style="margin-right:5px;" id="ressearch-field">
                <input id="ressearch-query" class="form-control" type="text" placeholder="Enter search keyword" aria-label="Search">
          </div>
          <div class="form-group" style="margin-right:5px;">
            <select id="ressearch-select" name="ressearch-select" class="selectpicker" data-width="fit">
              <option selected>in sequence</option>
              <option>TF genes</option>
              <option>p-value</option>
              <option>z-score</option>
            </select>
          </div>
          <div class="form-group">
            <button id="ressearch-btn" class="btn btn-primary btn-md" type="button">Search</button>
          </div>
          `);
        $('#ressearch-select').change(function(){
            var searchOpt = $("#ressearch-select option:selected").text();
            if(searchOpt == 'p-value' || searchOpt == 'z-score'){
                $('#math-compare-form').css("display","inline-block");
            }else{
                $('#math-compare-form').css("display","none");
            }
            if(searchOpt == 'TF genes'){
                $("#ressearch-field").html(`
                    <select id="ressearch-query" class="selectpicker" multiple data-live-search="true" data-actions-box="true"></select>
                `);
                var $rquery = $("#ressearch-query");
                $("#genes-dropdown > a").each(function(){
                    $rquery.append($('<option />').val($(this).text()).text($(this).text()));
                });
                $("#ressearch-query").selectpicker('refresh');
            }else{
                $("#ressearch-field").html(`
                    <input id="ressearch-query" class="form-control" type="text" placeholder="Enter search keyword" aria-label="Search" />
                `);
            }
        });
        $("#ressearch-query").keyup(function(event) {
            // Number 13 is the "Enter" key on the keyboard
            if (event.which === 13) {
                event.preventDefault();  // Cancel the default action, if needed
                // Trigger the button element with a click
                $("#ressearch-btn").click();
            }
        });
        $("#ressearch-btn").click(function() {
             // set search field for getrestbl
            var searchKey = $("#ressearch-query").val();
            var selected = $("#ressearch-select option:selected").text();
            var check = true;
            if(!searchKey || (selected == "TF genes" && $("#ressearch-query > option:selected").length==0)){
                check = false;
            }else if(isNaN(searchKey) && (selected == 'p-value' || selected == 'z-score')){
                check = false;
                alert("error: search value must be numeric"); // $('#search-error').html
            }
            if(check){
                if(selected == "p-value" || selected == "z-score"){
                    var compare = $("#math-compare-form option:selected").text();
                    addFilterElm(compare,searchKey,selected);
                }else if(selected == "TF genes"){
                    $("#ressearch-query > option:selected").each(function(){
                        addFilterElm("exact",$(this).text(),"TF_gene");
                    });
                }else{
                    addFilterElm(selected,searchKey,"sequence");
                }
                restbl.ajax.reload();
            }
        });
        $("#ressearch-select").selectpicker(); // to show the dropdownn options
        $("#math-compare-select").selectpicker();
    });
}

  /*
   * status_url:
   * parents: denotes the ids of the parent functions in the chain/pipeline
   */
function updateProgress(status_url,parents){
   $.getJSON(status_url,parents).done(function(data) {
       // update UI
       percent = parseInt(data['current'] * 100 / data['total']);

       $(".progress-bar").attr('aria-valuenow', percent).css('width',percent+"%").html(percent+"%"); // .attr('aria-valuenow', percent)
       $('#status').html(data['status']);

       if (data['state'] != 'PENDING' && data['state'] != 'PROGRESS') { // it's finish
           if ('result' in data) {
               //var res = data['result'];
               $('#status').hide();
               $(".progress").hide();
               displayResult(status_url);
               $('#csv-download').css('display','block').html("<a href=\"/files/" + data['taskid'] + "\" download=\"result.csv\">Download CSV</a>");
           }
           else if ('error' in data) {
               // found an error
               $('#status').html('Error: ' + data['error']);
           }else{
               // something unexpected happened
               $('#status').html('Result: ' + data['state']);
           }
       }
       else {
           // rerun in 2 seconds
           setTimeout(function() {
               updateProgress(status_url,parents);
           }, 2000);
       }
    });
 }

function inputProgress(status_url,parents){
    div = $('<div class="progress-bar" role="progressbar" style="width:0%"> \
                0% \
             </div>');
    // Progress bar -- Percentage -- Status message -- Result
    //$('#upload-form').html("<p>Processing your input...</p>");
    $('.progress').css('display','block').html(div);

    $('#status').css('display','block');
    updateProgress(status_url,JSON.parse(parents));
}

function getInputParam(status_url){
    var job_arr = status_url.split("/");
    var job_id = job_arr[job_arr.length-1];
    $.getJSON('/getinputparam/'+job_id,function(data) {
        var filterstr = "";

        /* String for desired output */
        if (data['filteropt'] == "1"){
            filterstr = "top " + data['filterval'] + " largest absolute z-score";
        }else{
            filterstr = "p-value < " + data['filterval'];
        }

        /* parse the returned list  */
        var $genesDropdown = $("#genes-dropdown");
        $.each(data['genes_selected'],function(idx,val){
            $genesDropdown.append("<a class=\"dropdown-item\" href=\"#\">"+val+"</a>");
        });
        $('.selectpicker').selectpicker('refresh');

        $('#inputpanel').prepend(
            '<p><b> File input:</b><br />' + data['filename'] + '</p>' +
            '<p><b> Desired output:</b><br />' + filterstr + '</p>' +
            '<p><b> Genome version:</b><br />' + data['chrver'] + '</p>'
        )
    });
}

$(function(){
    $('.navbar-nav .nav-link').removeClass('active');
    $('#nav-2').addClass('active');
});
