$( function() {
  
  var oldTop;
  var oldLeft;
  
  // the minimize code logic for UI Dialog boxes
  jQuery.fn.minimize = function() {
    
      var $this = $(this);
      $parent = $(this).parent().parent();
      
      if ($parent.hasClass("ui-dialog-minimized")) {
        var toRemove = $parent.parent();
        $("#gating_controls").dialog("option", "draggable", true);
        $parent.css("top", oldTop);
        $parent.css("left", oldLeft);
        $parent.appendTo("body");
        toRemove.remove();
      } else {
        oldTop = $parent.css("top");
        oldLeft = $parent.css("left");
        $("#gating_controls").dialog("option", "draggable", false);
        var $container = $("<div />", {
          style: "width: 320px; float: right;"
        });
        $parent.css("left", "").css("top", "");
        $parent.appendTo($container);
        $container.appendTo("#dialog-minimized-fixed-container");
      }
      
      $parent.toggleClass("ui-dialog-minimized");
      
  }
  
  // the function for generating a dialog box
  jQuery.fn.makeDialog = function() {
    
    // if the container for the bottom stuff doesn't exist, make it
    if ($("#dialog-minimized-fixed-container").length < 1) {
      var bot = $("<div />", {
        id: "dialog-minimized-fixed-container"
      });
      $("body").append(bot);
    }
    
    return this.each( function() {
      
      $(this).dialog({
        title: "Gating Method Controls",
        height: 500,
        minWidth: 250
      });
      
      var o = $(this).parent();
      
      // remove the 'close' button
      o.find("button").remove();
      
      // add 'minimize' button
      var btn = $("<button />", {
        class: "ui-button ui-widget ui-state-default ui-corner-all ui-button-icon-only ui-dialog-titlebar-minimize",
        role: "button",
        "aria-disabled": "false",
        title: "minimize"
      });
      
      btn
        .css("width", "22px")
        .css("height", "22px")
        .css("float", "right")
      ;
      
      var spn = $("<span />", {
        class: "ui-button-icon-primary ui-icon ui-icon-minus"
      });
      
      btn.append(spn);
      
      o.find(".ui-dialog-titlebar").append(btn);
      
    });
  }
  
  $("#gating_controls").makeDialog();
  $(".ui-dialog-titlebar-minimize").minimize();
  
  $(".ui-dialog-titlebar-minimize").click( jQuery.fn.minimize );
  
});
