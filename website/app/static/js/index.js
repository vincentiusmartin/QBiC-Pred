
function uploadPredResult(){
    var form = $('#predupload-form')[0];
    var fd = new FormData(form);
    $.ajax({
        type: 'POST',
        url: '/submitpredfile',
        data: fd, // pass the form data to Flask
        processData: false,  // tell jQuery not to process the data
        contentType: false,   // tell jQuery not to set contentType
        success: function(data, status, request) { // jsonified return from handle_upload in views.py
            status_url = request.getResponseHeader('Location');
            window.location.replace(status_url);
        },
        error: function(data,status) {
            alert(data["responseJSON"]["Message"]);
        }
    });
}

$(function() {
    $('#submit-pred').click(uploadPredResult);

    $('.navbar-nav .nav-link').removeClass('active');
    $('#nav-1').addClass('active');
});
