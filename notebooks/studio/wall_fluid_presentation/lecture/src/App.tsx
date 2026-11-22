/// <reference path="./marimo-studio.d.ts" />

import { Deck, Slide } from "@revealjs/react";
import { useEffect, useRef } from "react";
import type { RevealApi } from "reveal.js";
import "reveal.js/reveal.css";

const VIEW_TITLE = "Lecture";
const NOTEBOOK_TITLE = "Charged fluids near interfaces";

const keyboardCondition = (event: KeyboardEvent) =>
  !event.composedPath().some(
    (target) =>
      target instanceof Element && target.matches("marimo-cell, marimo-output"),
  );

/** Center slides again once projected notebook content has its final size. */
const useSettledLayout = () => {
  const deck = useRef<RevealApi | null>(null);
  useEffect(() => {
    const layout = () => deck.current?.layout();
    document.addEventListener("marimo-studio:idle", layout);
    return () => document.removeEventListener("marimo-studio:idle", layout);
  }, []);
  return deck;
};

export const App = () => {
  const deck = useSettledLayout();
  return (
    <Deck
      className="studio-deck"
      deckRef={deck}
      config={{
        controls: true,
        keyboardCondition,
        progress: true,
        scrollActivationWidth: 0,
        transition: "slide",
      }}
    >
      <Slide>
        <p className="deck-kicker">{VIEW_TITLE}</p>
        <h1>{NOTEBOOK_TITLE}</h1>
      </Slide>
      {[
        {
          "target": "title",
          "title": "Title",
          "showTitle": false,
        },
        {
          "target": "agenda",
          "title": "Agenda",
          "showTitle": false,
        },
        {
          "target": "problem",
          "title": "The problem",
          "showTitle": false,
        },
        {
          "target": "method",
          "title": "Method",
          "showTitle": false,
        },
        {
          "target": "stack",
          "title": "Software stack",
          "showTitle": false,
        },
        {
          "target": "fluids_table",
          "title": "Fluids",
          "showTitle": false,
        },
        {
          "target": "fig32_intro",
          "title": "Thesis Fig 3.2",
          "showTitle": false,
        },
        {
          "target": "fig32_plot",
          "title": "Wall profiles",
          "showTitle": false,
        },
        {
          "target": "fig32_stats",
          "title": "Peak summary",
          "showTitle": false,
        },
        {
          "target": "open_work",
          "title": "Open work",
          "showTitle": false,
        },
        {
          "target": "try_explorer",
          "title": "Try it",
          "showTitle": false,
        },
        {
          "target": "closing",
          "title": "Thanks",
          "showTitle": false,
        },
      ].map(({ target, title, showTitle }) => (
        <Slide key={target}>
          {showTitle ? <h2>{title}</h2> : null}
          <marimo-cell name={target} />
        </Slide>
      ))}
    </Deck>
  );
};
