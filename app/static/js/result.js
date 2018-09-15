
function displayResult(status_url){
    var keys = [];
    // show result
    var job_arr = status_url.split("/");
    var job_id = job_arr[job_arr.length-1];
    var cols;
    $.getJSON('/getrescol/'+job_id,function(cols) {
        searchKey = "";
        searchOpt = "All";

        $('#tablewrap').css('display','block');
        var restbl = $('#restbl').DataTable({
            dom: 'lr<"#ressearch">tip',
            searching: false,
            serverSide: true,
            orderMulti: false,
            ajax: {
              url: '/getrestbl/'+job_id,
              data:function(d){
                 d.csKey = searchKey; // custom search
                 d.csOpt = searchOpt;
              }
            },
            columns: [
               { title: cols[0] },
               { title: cols[1] },
               { title: cols[2] },
               { title: cols[3] },
               { title: cols[4] },
               { title: cols[5] },
            ]
        });
        // This is used to replace the search bar from DataTable so we can have
        // our custom. Search bar itself is inactivated.
        // https://stackoverflow.com/questions/43454575/call-datatable-ajax-call-on-custom-button
        $("#ressearch").addClass("form-inline mb-3").html(`
          <div class="form-grup" style="margin-right:5px;">
            <input id="ressearch-text" class="form-control" type="text" placeholder="Enter search keyword" aria-label="Search">
          </div>
          <div class="form-grup" style="margin-right:5px;">
            <select id="ressearch-select" name="ressearch-select" class="selectpicker" data-width="fit">
              <option selected>All</option>
              <option>Opt-233</option>
            </select>
          </div>
          <div class="form-grup">
           <button id="ressearch-btn" class="btn btn-primary btn-md" type="button">Search</button>
          </div>
          `);
         $("#ressearch-btn").click(function() {
            searchKey = $("#ressearch-text").val();
            searchOpt = $("#ressearch-select option:selected").text();
            restbl.ajax.reload();
         });
         $("#ressearch-select").selectpicker(); // to show the dropdownn options
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
                   $('#csv-download').css('display','block').html("<a href=\"/files/" + data['csvlink'] + "\" download=\"result.csv\">Download CSV</a>");
               }
               else {
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
        $.each(data['pbmselected'],function(idx,val){
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
