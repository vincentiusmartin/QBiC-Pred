
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
                 var res = data['result'];
                 var tbl = '<table class="table table-bordered table-hover table-striped" id="restbl">';
                 var keys = [];
                 tbl += '<thead><tr>';
                 for (var i = 0; i < Object.keys(res[0]).length; i++) { // zero index is the list of colname
                     keys.push(res[0][i]);
                     tbl += '<th>' + res[0][i] + '</th>';
                 }
                 tbl += '</tr></thead>';
                 tbl += '<tbody>';
                 for (var i = 1; i < res.length; i++) {
                     tbl += '<tr>';
                     for (var j = 0; j < keys.length; j++){
                         tbl += '<td>' + res[i][keys[j]] + '</td>';
                     }
                     tbl += '</tr>';
                 }
                 tbl += '</tbody>';
                 tbl += '</table>';
                 // show result
                 $('#status').html('Result: ' + tbl);
                 var job_arr = status_url.split("/");
                 $('#restbl').DataTable({
                   searching: false,
                   info: false
                  });
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
