$( function() {
  $("#R_input").keydown( function(evt) {
    if (evt.keyCode == 13) { // enter
      $("#R_send").trigger("click");
    }
  });
})
