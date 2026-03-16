// fmin.xyz Touying presentation theme for Quarto
#import "@preview/touying:0.6.1": *
#import "@preview/fontawesome:0.5.0": *

// Main accent colour
#let main_colour = rgb(51, 51, 179)

// Override Quarto default callout with our theme style
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black, body_background_color: white) = {
  block(
    width: 100%,
    fill: rgb("#f0f2f8"),
    stroke: (left: 2pt + main_colour.lighten(20%)),
    inset: (x: 6pt, y: 5pt),
    radius: (right: 2pt),
  )[
    #if title != none and title != [] [
      #text(weight: "bold", size: 9pt, title)
      #v(1pt)
    ]
    #body
  ]
}

// Common math operators
#let argmin = math.op("argmin", limits: true)
#let argmax = math.op("argmax", limits: true)
#let sign = math.op("sign")
#let dom = math.op("dom")
#let conv = math.op[#strong[conv]]
#let ri = math.op[#strong[ri]]
#let Span = math.op("span")
#let prox = math.op("prox")

// Helper: image from files/ directory
#let img(path, ..args) = image("../files/" + path, ..args)

// Helper: wrap a list/enum to make items appear one by one
#let incremental(body) = {
  if type(body) != content { return body }
  let children = if body.has("children") { body.children } else { (body,) }
  let items = children.filter(c => type(c) == content and (c.func() == list.item or c.func() == enum.item))
  if items.len() > 1 { items.join(pause) } else { body }
}

// ============================================================
// Footer
// ============================================================

#let _footer-overlay = {
  context {
    let ratio = utils.slide-counter.get().first() / utils.slide-counter.final().first()
    place(bottom + left, rect(width: ratio * 100%, height: 0.2pt, fill: rgb(255, 102, 0)))
  }
  place(bottom + left,
    link("https://fmin.xyz", image("../files/logo.svg", height: 0.30cm)))
  place(bottom, dy: -0.05cm, box(width: 100%, inset: (left: 1.3cm, right: 0.08cm), {
    set text(size: 6pt, font: "CMU Sans Serif", weight: "regular", style: "normal")
    context {
      let cp = here().page()
      let hs = query(heading).filter(h => h.location().page() <= cp and h.level == 1)
      if hs.len() > 0 {
        link(hs.last().location(), text(size: 5pt, font: "CMU Sans Serif", weight: "regular", hs.last().body))
      }
    }
    h(1fr)
    link("https://mipt25.fmin.xyz", text(size: 5pt, font: "Font Awesome 6 Free", weight: "regular", "\u{f3a5}"))
    h(0.08cm)
    link("https://github.com/MerkulovDaniil/mipt25", text(size: 5pt, font: "Font Awesome 6 Brands", "\u{f09b}"))
    h(0.1cm)
    link("https://t.me/fminxyz", text(size: 5pt, font: "Font Awesome 6 Brands", "\u{f2c6}"))
    h(0.3cm)
    box(width: 0.4cm, align(right, context utils.slide-counter.display()))
  }))
}

// ============================================================
// Slide types
// ============================================================

#let slide(
  config: (:),
  repeat: auto,
  setting: body => body,
  composer: auto,
  ..bodies,
) = touying-slide-wrapper(self => {
  let header(self) = {
    set align(left + horizon)
    if utils.display-current-heading(level: 2) != [] {
      block(
        width: 100%,
        height: 0.85cm,
        inset: (x: 0pt),
        context {
          let heading-content = utils.display-current-heading(level: 2)
          let base-size = 12pt
          let avail = 160mm - 2 * 0.214cm
          let m = measure(text(size: base-size, weight: "bold", style: "italic", font: "CMU Sans Serif", heading-content))
          let sz = if m.width > avail { base-size * (avail / m.width) } else { base-size }
          text(
            size: sz,
            weight: "bold",
            style: "italic",
            fill: main_colour,
            font: "CMU Sans Serif",
            heading-content,
          )
        },
      )
    }
  }
  let self = utils.merge-dicts(
    self,
    config-page(header: header),
  )
  touying-slide(
    self: self,
    config: config,
    repeat: repeat,
    setting: body => setting(block(clip: true, height: 90mm - 0.85cm - 0.45cm - 0.85cm, width: 100%, body)),
    composer: composer,
    ..bodies,
  )
})

#let title-slide(
  config: (:),
  bg-image: none,
) = touying-slide-wrapper(self => {
  let bg = if bg-image != none { bg-image } else { self.store.bg-image }
  let self = utils.merge-dicts(
    self,
    config-page(
      header: self => [],
    ),
  )
  touying-slide(self: self, config: config, {
    place(center + horizon, dx: -0.214cm, dy: -0.2cm, image(bg, width: 170mm, height: 96mm))
    place(
      center + horizon,
      block(
        width: 90mm,
        fill: rgb(255, 255, 255, 70%),
        radius: 10pt,
        inset: (x: 10pt, y: 8pt),
      )[
        #set align(center)
        #set text(font: "CMU Sans Serif")
        #text(weight: "bold", size: 11pt, self.info.title)
        #v(0.35cm)
        #text(weight: "bold", size: 11pt, self.info.author)
        #v(0.35cm)
        #text(weight: "bold", size: 8pt, self.info.institution)
      ],
    )
  })
})

#let new-section-slide(config: (:), body) = touying-slide-wrapper(self => {
  let self = utils.merge-dicts(
    self,
    config-page(
      header: self => [],
    ),
  )
  touying-slide(
    self: self,
    config: config,
    align(
      center + horizon,
      text(size: 20pt, weight: "bold", fill: main_colour, font: "CMU Sans Serif", utils.display-current-heading(level: 1)),
    ),
  )
})

// ============================================================
// Helper: theorem callout block
// ============================================================

#let theorem-block(title: none, body) = {
  block(
    width: 100%,
    fill: rgb("#f0f2f8"),
    stroke: (left: 2pt + main_colour.lighten(20%)),
    inset: (x: 6pt, y: 5pt),
    radius: (right: 2pt),
  )[
    #if title != none [
      #text(weight: "bold", title)
      #v(1pt)
    ]
    #body
  ]
}

// ============================================================
// Main theme function
// ============================================================

#let fmin-theme(
  aspect-ratio: "16-9",
  handout: false,
  bg-image: "../files/back19.jpeg",
  ..args,
  body,
) = {
  show: touying-slides.with(
    config-page(
      width: 160mm,
      height: 90mm,
      margin: (x: 0.214cm, top: 0.85cm, bottom: 0.45cm),
      footer: none,
      footer-descent: 0em,
      header-ascent: 0em,
      foreground: _footer-overlay,
    ),
    config-common(
      handout: handout,
      frozen-counters: (counter(figure),),
      slide-fn: slide,
      new-section-slide-fn: new-section-slide,
      slide-level: 2,
      zero-margin-header: false,
      zero-margin-footer: false,
    ),
    config-methods(
      init: (self: none, body) => {
        set text(
          font: "CMU Sans Serif",
          size: 9pt,
          lang: "ru",
        )
        set par(leading: 0.4em, spacing: 0.65em)
        show math.equation: set text(font: "New Computer Modern Math")
        set math.equation(numbering: "(1)", supplement: none)
        show math.equation.where(block: true): it => {
          if it.has("label") or it.numbering == none { it }
          else {
            counter(math.equation).update(v => v - 1)
            math.equation(block: true, numbering: none, it.body)
          }
        }
        set list(marker: (sym.circle.filled, sym.circle.filled, sym.circle.filled))
        set enum(numbering: "1.")
        set table(stroke: none, inset: 6pt, align: center + horizon)
        set figure(supplement: [Рис.])
        set figure.caption(separator: [: ])
        show figure: set block(breakable: false)
        body
      },
    ),
    config-colors(
      primary: main_colour,
      neutral-lightest: rgb("#ffffff"),
      neutral-darkest: rgb("#000000"),
    ),
    config-store(
      bg-image: bg-image,
    ),
    config-info(..args.named()),
    ..args,
  )
  body
}
