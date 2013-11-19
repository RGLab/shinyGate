$( function() {
  $("#R_input").keydown( function(evt) {
    if (evt.keyCode == 13) {
      $("#R_send").trigger("click");
    }
  });
})
