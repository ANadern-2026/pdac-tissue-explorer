# Minimal placeholder for the retired Connect Cloud address: forwards visitors to the app's permanent home.
# Published by deploy/connect-cloud/retire.R in place of the full app (no data, very little memory).
library(shiny)
NEW_URL <- "https://pdac.biointeract.space/"

ui <- fluidPage(
  tags$head(
    tags$title("PDAC Tissue Explorer has moved"),
    tags$meta(`http-equiv` = "refresh", content = paste0("3; url=", NEW_URL)),
    tags$script(HTML(sprintf("setTimeout(function () { try { window.top.location.replace('%s'); } catch (e) { window.location.replace('%s'); } }, 1500);", NEW_URL, NEW_URL))),
    tags$style(HTML("body { font-family: system-ui, -apple-system, 'Segoe UI', Roboto, Arial, sans-serif; background: #fcfcfb; color: #0b0b0b; }
                     .box { max-width: 560px; margin: 12vh auto; padding: 2rem; border: 1px solid #e9e8e4; border-radius: 10px; background: #fff; }
                     a.btn-go { display: inline-block; margin-top: 0.8rem; padding: 0.55rem 1.1rem; background: #2a78d6; color: #fff; border-radius: 6px; text-decoration: none; }"))),
  div(class = "box",
      h2("PDAC Tissue Explorer has moved"),
      p("The app now runs at its permanent address. You will be forwarded automatically."),
      a(class = "btn-go", href = NEW_URL, target = "_top", "Open pdac.biointeract.space"),
      p(style = "color: #52514e; font-size: 0.85rem; margin-top: 1.2rem;", "Please update your bookmarks.")))

shinyApp(ui, function(input, output, session) {})
