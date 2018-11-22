
function testGetTbl(ajax_url){
    // table eneds to be initialized in this format otherwise it wouldn't work
    var tbl =
    `<table class="table table-bordered table-hover table-striped" id="restbl">
    <thead>
      <tr>
       <th></th>
       <th></th>
       <th></th>
       <th></th>
       <th></th>
       <th></th>
      </tr>
     </thead>
     <tbody>
      <tr>
       <td></td>
       <td></td>
       <td></td>
       <td></td>
       <td></td>
       <td></td>
      </tr>
    </tbody>
    </table>`;
    var cols;
    $.getJSON('/gettestcol',function(cols) {
        $('#testingmain').html('Result: ' + tbl);
        $('#restbl').DataTable({
            serverSide: true,
            ajax: {
              url: ajax_url,
            },
            columns: [
               { title: cols[0] },
               { title: cols[1] },
               { title: cols[2]  },
               { title: cols[3]  },
               { title: cols[4] },
               { title: cols[5]  },
            ]
        });
    });

    /*$.getJSON(url).done(function(datatbl) {
        var tbl = '<table class="table table-bordered table-hover table-striped" id="restbl"></table>';
        $('#testingmain').html('Result: ' + tbl);
        $('#restbl').DataTable({
         data: datatbl["data"],
         columns: [
            { data: datatbl["cols"][0], title: datatbl["cols"][0]},
            { data: datatbl["cols"][1], title: datatbl["cols"][1]},
            { data: datatbl["cols"][2], title: datatbl["cols"][2]  },
            { data: datatbl["cols"][3], title: datatbl["cols"][3]  },
            { data: datatbl["cols"][4], title: datatbl["cols"][4] },
            { data: datatbl["cols"][5], title: datatbl["cols"][5]  },
         ]
        });
    });*/
}
