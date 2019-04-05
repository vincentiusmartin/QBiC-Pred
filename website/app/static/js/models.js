

function fillModelsTbl(){
    // TypeError: e[i] is undefined: means that column names are not there
    $.getJSON('/getmodeltbl',function(tbl) {
        if(tbl["data"].length > 0){
            var colNames = []
            for(var i = 0; i < tbl["cols"].length; i++){
                // need to specify both title and data for colnames
                colNames.push({"title":tbl["cols"][i], "data":tbl["cols"][i]});
            }
            var restbl = $('#models-tbl').DataTable({
                data: tbl["data"],
                columns: colNames
            });
        }
    });
    //$('#models-tbl').html(fill);
    //$('#models-tbl').DataTable();
}

$(function(){
    //displayModelsTbl();
    fillModelsTbl();

    $('.navbar-nav .nav-link').removeClass('active');
    $('#nav-4').addClass('active');
});
