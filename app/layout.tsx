import type { Metadata } from 'next'
import './globals.css'

export const metadata: Metadata = {
  title: 'rRNA Operon Evolution Analysis',
  description: 'Neutral model analysis for rRNA operon evolution',
}

export default function RootLayout({
  children,
}: {
  children: React.ReactNode
}) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  )
}