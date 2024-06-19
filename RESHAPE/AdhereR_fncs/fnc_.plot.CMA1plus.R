############################################################################################
#
# Plot CMAs greater or equal to 1 (smart enough to know what to do with specific info)
#
############################################################################################
.plot.CMA1plus <- function(cma,                                   # the CMA1 (or derived) object
                           patients.to.plot=NULL,                 # list of patient IDs to plot or NULL for all
                           duration=NA,                           # duration to plot in days (if missing, determined from the data)
                           align.all.patients=FALSE, align.first.event.at.zero=FALSE, # should all patients be aligned? and, if so, place the first event as the horizintal 0?
                           show.period=c("dates","days")[2],      # draw vertical bars at regular interval as dates or days?
                           period.in.days=90,                     # the interval (in days) at which to draw veritcal lines
                           show.legend=TRUE, legend.x="right", legend.y="bottom", legend.bkg.opacity=0.5, legend.cex=0.75, legend.cex.title=1.0, # legend params and position
                           cex=1.0, cex.axis=0.75, cex.lab=1.0,   # various graphical params
                           show.cma=TRUE,                         # show the CMA type
                           xlab=c("dates"="Date", "days"="Days"), # Vector of x labels to show for the two types of periods, or a single value for both, or NULL for nothing
                           ylab=c("withoutCMA"="patient", "withCMA"="patient (& CMA)"), # Vector of y labels to show without and with CMA estimates, or a single value for both, or NULL ofr nonthing
                           title=c("aligned"="Event patterns (all patients aligned)", "notaligned"="Event patterns"), # Vector of titles to show for and without alignment, or a single value for both, or NULL for nothing
                           col.cats=rainbow,                      # single color or a function mapping the categories to colors
                           unspecified.category.label="drug",     # the label of the unspecified category of medication
                           medication.groups.to.plot=NULL,        # the names of the medication groups to plot (by default, all)
                           medication.groups.separator.show=TRUE, medication.groups.separator.lty="solid", medication.groups.separator.lwd=2, medication.groups.separator.color="blue", # group medication events by patient?
                           medication.groups.allother.label="*",  # the label to use for the __ALL_OTHERS__ medication class (defaults to *)
                           lty.event="solid", lwd.event=2, pch.start.event=15, pch.end.event=16, # event style
                           show.event.intervals=TRUE,             # show the actual prescription intervals
                           plot.events.vertically.displaced=TRUE, # display the events on different lines (vertical displacement) or not (defaults to TRUE)?
                           col.na="lightgray",                    # color for missing data
                           print.CMA=TRUE, CMA.cex=0.50,          # print CMA next to the participant's ID?
                           plot.CMA=TRUE,                   # plot the CMA next to the participant ID?
                           CMA.plot.ratio=0.10,             # the proportion of the total horizontal plot to be taken by the CMA plot
                           CMA.plot.col="lightgreen", CMA.plot.border="darkgreen", CMA.plot.bkg="aquamarine", CMA.plot.text=CMA.plot.border, # attributes of the CMA plot
                           highlight.followup.window=TRUE, followup.window.col="green",
                           highlight.observation.window=TRUE, observation.window.col="yellow", observation.window.density=35, observation.window.angle=-30, observation.window.opacity=0.3,
                           show.real.obs.window.start=TRUE, real.obs.window.density=35, real.obs.window.angle=30, # for some CMAs, the real observation window starts at a different date
                           print.dose=FALSE, cex.dose=0.75, print.dose.outline.col="white", print.dose.centered=FALSE, # print daily dose
                           plot.dose=FALSE, lwd.event.max.dose=8, plot.dose.lwd.across.medication.classes=FALSE, # draw daily dose as line width
                           alternating.bands.cols=c("white", "gray95"), # the colors of the alternating vertical bands across patients (NULL=don't draw any; can be >= 1 color)
                           bw.plot=FALSE,                         # if TRUE, override all user-given colors and replace them with a scheme suitable for grayscale plotting
                           rotate.text=-60,                       # some text (e.g., axis labels) may be rotated by this much degrees
                           force.draw.text=FALSE,                 # if true, always draw text even if too big or too small
                           min.plot.size.in.characters.horiz=0, min.plot.size.in.characters.vert=0, # the minimum plot size (in characters: horizontally, for the whole duration, vertically, per event)
                           suppress.warnings=FALSE,               # suppress warnings?
                           max.patients.to.plot=100,              # maximum number of patients to plot
                           export.formats=NULL,                   # the formats to export the figure to (by default, none); can be any subset of "svg" (just SVG file), "html" (SVG + HTML + CSS + JavaScript all embedded within the HTML document), "jpg", "png", "webp", "ps" and "pdf"
                           export.formats.fileprefix="AdhereR-plot", # the file name prefix for the exported formats
                           export.formats.height=NA, export.formats.width=NA, # desired dimensions (in pixels) for the exported figure (defaults to sane values)
                           export.formats.save.svg.placeholder=TRUE,
                           export.formats.svg.placeholder.type=c("jpg", "png", "webp")[2],
                           export.formats.svg.placeholder.embed=FALSE, # save a placeholder for the SVG image?
                           export.formats.html.template=NULL, export.formats.html.javascript=NULL, export.formats.html.css=NULL, # HTML, JavaScript and CSS templates for exporting HTML+SVG
                           export.formats.directory=NA,           # if exporting, which directory to export to (if not give, creates files in the temporary directory)
                           generate.R.plot=TRUE,                  # generate standard (base R) plot for plotting within R?
                           do.not.draw.plot=FALSE,                # if TRUE, don't draw the actual plot, but only the legend (if required)
                           ...
)
{
  .plot.CMAs(cma,
             patients.to.plot=patients.to.plot,
             duration=duration,
             align.all.patients=align.all.patients,
             align.first.event.at.zero=align.first.event.at.zero,
             show.period=show.period,
             period.in.days=period.in.days,
             show.legend=show.legend,
             legend.x=legend.x,
             legend.y=legend.y,
             legend.bkg.opacity=legend.bkg.opacity,
             legend.cex=legend.cex,
             legend.cex.title=legend.cex.title,
             cex=cex,
             cex.axis=cex.axis,
             cex.lab=cex.lab,
             show.cma=show.cma,
             xlab=xlab,
             ylab=ylab,
             title=title,
             col.cats=col.cats,
             unspecified.category.label=unspecified.category.label,
             medication.groups.to.plot=medication.groups.to.plot,
             medication.groups.separator.show=medication.groups.separator.show,
             medication.groups.separator.lty=medication.groups.separator.lty,
             medication.groups.separator.lwd=medication.groups.separator.lwd,
             medication.groups.separator.color=medication.groups.separator.color,
             medication.groups.allother.label=medication.groups.allother.label,
             lty.event=lty.event,
             lwd.event=lwd.event,
             show.event.intervals=show.event.intervals,
             plot.events.vertically.displaced=plot.events.vertically.displaced,
             pch.start.event=pch.start.event,
             pch.end.event=pch.end.event,
             print.dose=print.dose,
             cex.dose=cex.dose,
             print.dose.outline.col=print.dose.outline.col,
             print.dose.centered=print.dose.centered,
             plot.dose=plot.dose,
             lwd.event.max.dose=lwd.event.max.dose,
             plot.dose.lwd.across.medication.classes=plot.dose.lwd.across.medication.classes,
             col.na=col.na,
             print.CMA=print.CMA,
             CMA.cex=CMA.cex,
             plot.CMA=plot.CMA,
             CMA.plot.ratio=CMA.plot.ratio,
             CMA.plot.col=CMA.plot.col,
             CMA.plot.border=CMA.plot.border,
             CMA.plot.bkg=CMA.plot.bkg,
             CMA.plot.text=CMA.plot.text,
             plot.partial.CMAs.as=FALSE, # CMA1+ do not have partial CMAs
             highlight.followup.window=highlight.followup.window,
             followup.window.col=followup.window.col,
             highlight.observation.window=highlight.observation.window,
             observation.window.col=observation.window.col,
             observation.window.density=observation.window.density,
             observation.window.angle=observation.window.angle,
             observation.window.opacity=observation.window.opacity,
             show.real.obs.window.start=show.real.obs.window.start,
             real.obs.window.density=real.obs.window.density,
             real.obs.window.angle=real.obs.window.angle,
             alternating.bands.cols=alternating.bands.cols,
             bw.plot=bw.plot,
             rotate.text=rotate.text,
             force.draw.text=force.draw.text,
             min.plot.size.in.characters.horiz=min.plot.size.in.characters.horiz,
             min.plot.size.in.characters.vert=min.plot.size.in.characters.vert,
             max.patients.to.plot=max.patients.to.plot,
             export.formats=export.formats,
             export.formats.fileprefix=export.formats.fileprefix,
             export.formats.height=export.formats.height,
             export.formats.width=export.formats.width,
             export.formats.save.svg.placeholder=export.formats.save.svg.placeholder,
             export.formats.svg.placeholder.type=export.formats.svg.placeholder.type,
             export.formats.svg.placeholder.embed=export.formats.svg.placeholder.embed,
             export.formats.html.template=export.formats.html.template,
             export.formats.html.javascript=export.formats.html.javascript,
             export.formats.html.css=export.formats.html.css,
             export.formats.directory=export.formats.directory,
             generate.R.plot=generate.R.plot,
             do.not.draw.plot=do.not.draw.plot,
             suppress.warnings=suppress.warnings);
}