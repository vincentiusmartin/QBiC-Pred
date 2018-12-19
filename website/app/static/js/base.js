
function getRecentJobs(){
    var $dropdown = $("#dropdown-header");
    $.ajax({
        type: "GET",
        url: '/recent',
        dataType: 'json',
        success: function(data){
            $.each(data, function(k, v) {
                $dropdown.append("<a class=\"dropdown-item\" href=\"" + v[1] + "\">" + v[0] + "</a>");
            });
        },
        error: function() {
          alert("error getting recent jobs");
       }
    });
}


$(function() {
    getRecentJobs();
});
